# RRTF1-dependent_touch_transcriptome
Please cite: Park, S., Finlayson, S.A., C, Li., 2026. RRTF1 promotes touch-responses in Arabidopsis shoots independent of jasmonic acid. [BioRxiv](doi: https://doi.org/10.64898/2026.03.08.710212).

# RNA-seq analysis of Arabidopsis shoots.
1. Trimmomatic_code.sh runs Trimmomatic to trim the adapter sequences and low-quality reads.
2. HISAT2_code_at.sh runs HISAT2 to align the trimmed paired-end and unpaired reads to the Arabidopsis thaliana reference genome, coverts the output to sorted BAM files, and filters for reads alignments with a mapping quality (MAPQ) score of 30 or higher.
3. HISAT2_build_index_AT runs HISAT2 to build the Arabidopsis thaliana reference genome.
4. featureCounts_AT.sh runs featureCounts to quantify the number of reads mapped to each gene using the TAIR10 GTF annotation and the quality-filtered BAM files, generating a count matrix.

# Differential gene expression in R.
- For a time-course analysis of Col-0 was performed at 0, 5, 10, and 30 min post-touch, BAM files obtained from sequencing lanes 1 and 2 were merged before conducting differential gene expression analysis.
- DESeq_touch_response.R uses DESeq2 to analyze differential gene expression for a time-course analysis of Col-0 in response to touch.
- DESeq_rrtf1_touch.R uses DESeq2 to analyze differential gene expression between Col-0 and the rrtf1 mutant at 0 and 10 min post-touch. This script also integrates pooled RNA-seq RNA-seq data and performs batch correction.
- edgeR_rrtf1_touch.R uses edgeR to analyze differential gene expression using the glmTreat function.

# Figure generation
- Touch_figure.R generates figures.
- Touch_figure_supp.R generates supplementary figures.
