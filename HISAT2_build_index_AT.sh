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
hisat2-build -p 28 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz at_index_55

###splicesites.txt ###
###Python load
module load GCCcore/12.2.0
module load Python/3.10.8
hisat2_extract_splice_sites.py Arabidopsis_thaliana.TAIR10.55.gtf.gz > at_splicesites.txt