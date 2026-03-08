#!/bin/bash
#SBATCH --export=NONE
#SBATCH --job-name=job_name
#SBATCH --time=0-06:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=56G
#SBATCH --output=stdout.%x.%j
#SBATCH --error=stderr.%x.%j

module load GCC/11.2.0  OpenMPI/4.1.1 HISAT2/2.2.1
module load SAMtools/1.14
################################### VARIABLES ##################################

##Parameters##
threads=20
#read group information
id='af_amp'
library='sra'
platform='ILLUMINA'
genome_index_prefix='/scratch/user/psk/thigmo_at/genome_index/at_index_55'

for i in $(ls  *_1P.fastq.gz | sed 's/_trim_1P.fastq.gz//g' | sort -u)
do 
output_bam="${i}_pe.bam"
hisat2 --dta --known-splicesite-infile '/scratch/user/psk/thigmo_at/genome_index/at_splicesites.txt' -p $threads --rg-id "$id" --rg "LB:$library" --rg "SM:${i}" --rg "PL:$platform" -x $genome_index_prefix -q -1 "${i}_trim_1P.fastq.gz" -2 "${i}_trim_2P.fastq.gz" -U  ${i}_trim_1U.fastq.gz,${i}_trim_2U.fastq.gz -S out.sam
samtools view -bS out.sam | samtools sort -m 2G -@ $threads - -T $i -o $output_bam
samtools view -@ 20 -q 30 -b $output_bam > ${i}_q30.bam
rm out.sam
done