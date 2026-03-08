# RRTF1-dependent_touch_transcriptome
Please cite: Park, S., Finlayson, S.A., C, Li., 2026. RRTF1 promotes touch-responses in Arabidopsis shoots independent of jasmonic acid. BioRxiv.

# RNA-seq analysis of Arabidopsis shoots.
1. Trimmomatic_code.sh runs Trimmomatic to trim the adapter sequences and low-quality reads.
2. HISAT2_code_at.sh runs HISAT2 to align the trimmed paired-end and unpaired reads to the Arabidopsis thaliana reference genome, coverts the output to sorted BAM files, and filters for reads alignments with a mapping quality (MAPQ) score of 30 or higher.
3. HISAT2_build_index_AT runs HISAT2 to build the Arabidopsis thaliana reference genome.
4. featureCounts_AT.sh runs featureCounts to quantify the number of reads mapped to each gene using the TAIR10 GTF annotation and the quality-filtered BAM files, generating a count matrix.

# Differential gene expression in R.
1. 
