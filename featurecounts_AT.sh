#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=trimmomatic      # job name
#SBATCH --time=00-06:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=2           # CPUs (threads) per command
#SBATCH --mem=2G                    # total memory per node
#SBATCH --output=stdout.%j          # save stdout to file
#SBATCH --error=stderr.%j           # save stderr to file

#run together L1 and L2
module load GCC/11.2.0 Subread/2.0.3
featureCounts -p -t exon -g gene_id -a '/scratch/user/psk/thigmo_at/genome_index/Arabidopsis_thaliana.TAIR10.55.gtf' -o counts_at_q30.txt *_q30.bam