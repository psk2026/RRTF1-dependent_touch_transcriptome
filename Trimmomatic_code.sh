#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=trimmomatic      # job name
#SBATCH --time=00-06:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=2           # CPUs (threads) per command
#SBATCH --mem=2G                    # total memory per node
#SBATCH --output=stdout.%j          # save stdout to file
#SBATCH --error=stderr.%j           # save stderr to file

module load Trimmomatic/0.39-Java-11
################################### VARIABLES ##################################
# TODO Edit these variables as needed:
########## INPUTS ##########
threads=2
quality_format="-phred33"
adapter_file='TruSeq3-PE-2.fa'

## COMMANDS ##
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $threads $quality_format $pe1_1 $pe1_2 -baseout ${prefix}_trim.fastq.gz ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/$adapter_file:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36