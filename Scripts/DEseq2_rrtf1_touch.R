#####DESeq_for_Merged_Bam_file
setwd("~/../Thigmo_exp/2025_RNA_seq/rrtf1_RNAseq_033125/")

###First, do check raw count files.####
counts <- read.table("counts_rrtf_at_q30.txt",header=TRUE,row.names = 1)[,-c(1:5)]
at_metadata <- read.table("at_metaData.csv",header=TRUE,row.names = 1,sep=',')
at_metadata$time <- factor(at_metadata$time)
at_metadata$geno <- factor(at_metadata$geno)
# Check reference factor levels, if not 'con' for treatment and '1' for time relevel
print(at_metadata$time)
at_metadata$geno <- relevel(at_metadata$geno, ref = 'rrtf')
# next three steps should be in exact order 
sample_order <- order(at_metadata$geno)
at_metadata <- at_metadata[sample_order, , drop=F] #drop=F ensures that even if you're subsetting a single row, the result will still be a data frame,not a list.
counts <- counts[, rownames(at_metadata)]
colnames(counts)
rownames(at_metadata)
colnames(counts) == rownames(at_metadata)

library(dplyr)
library(ggplot2)
library(DESeq2)
# Use multiple cores for DEseq2
library("BiocParallel")
register(MulticoreParam(8))
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = at_metadata,
                              design = ~ time + geno)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = at_metadata,
                              design = ~ merge)
keep <- rowSums(counts(dds) > 10) >= 3
dds_keep <- dds[keep,]
#Normalized counts
dds_keep <- DESeq(dds_keep,betaPrior=FALSE, parallel=TRUE)
resultsNames(dds_keep)
# Check sample outlier using Cooks distance
summary(results(dds_keep, alpha=0.05))
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_keep)[["cooks"]]), range=0, las=2)
#Check within condition variability using PCA
vsd <- vst(dds_keep)
plotPCA(vsd, intgroup=c("time","geno"))
y<-plotPCA(vsd, intgroup=c("time","geno"))
y + geom_text(aes(label=substr(name, start = 1, stop = 8)),vjust=2,check_overlap = T,size = 2)
cbind(y$data$treatment,y$data$time,y$data$PC1,y$data$PC2)

at_metadata <- read.table("at_metaData.csv",header=TRUE,row.names = 1,sep=',')
at_metadata$merge <- factor(at_metadata$merge)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = at_metadata,
                              design = ~ merge)
keep <- rowSums(counts(dds) > 10) >= 3
dds_keep <- dds[keep,]
#Normalized counts
dds_keep <- DESeq(dds_keep,betaPrior=FALSE, parallel=TRUE)
resultsNames(dds_keep)
# Check sample outlier using Cooks distance
summary(results(dds_keep, alpha=0.05))
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_keep)[["cooks"]]), range=0, las=2)
#Check within condition variability using PCA
vsd <- vst(dds_keep)
plotPCA(vsd, intgroup=c("merge"))
y<-plotPCA(vsd, intgroup=c("merge"))
y + geom_text(aes(label=substr(name, start = 1, stop = 8)),vjust=2,check_overlap = T,size = 2)
cbind(y$data$treatment,y$data$time,y$data$PC1,y$data$PC2)

###PCA for rrtf1
library(DESeq2)    
library(ggplot2)    
library(dplyr)     
library(ggrepel)   
library(ggnewscale) 

pca_df <- plotPCA(vsd,                     
                  intgroup=c("merge"),
                  returnData = TRUE)
percent <- round(100 * attr(pca_df, "percentVar"))
pca_df$time <- as.factor(pca_df$time)
ggplot(pca_df, aes(x = PC1, y = PC2,
                   colour = group)) +
  scale_colour_manual(values = c(
    "rrtf0" = "red",    # vivid red
    "wt0" = "green",   # bright blue
    "rrtf10" = "orange",   # strong orange
    "wt10" = "blue"  # bright purple
  )) +
  geom_hline(yintercept = 0, linewidth = .3) +  # 
  geom_vline(xintercept = 0, linewidth = .3) +  # 
  geom_point(size = 4, stroke = 1) +
  labs(x = sprintf("PCA1: %s%% variance", percent[1]),
       y = sprintf("PCA2: %s%% variance", percent[2]),
       colour = "Group") +    # 
  theme_bw(base_size = 18) +
  theme(legend.position = "right",
        panel.grid = element_blank())


stat_ellipse(aes(group = group), linewidth = .7,
             linetype = "solid", alpha = .15) +
  
  
  ## PCA plot
  plotPCA(d,cls = 'day')

t <- data.frame(PC1 = y$data$PC1, PC2 = y$data$PC2, Group = merge)
library(ggplot2)
ggplot(pca_data_plot, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(variance_explained[1], 0), "% variance")) +
  ylab(paste0("PC2: ", round(variance_explained[2], 0), "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal() # or other themes

### Using all samples / get the model matrix
mod_mat <- model.matrix(design(dds_keep), colData(dds_keep))

# Define coefficient vectors for each condition
wt0 <- colMeans(mod_mat[dds_keep$merge == "wt0", ])
wt10 <- colMeans(mod_mat[dds_keep$merge == "wt10", ])
rrtf0 <- colMeans(mod_mat[dds_keep$merge == "rrtf0", ])
rrtf10 <- colMeans(mod_mat[dds_keep$merge == "rrtf10", ])

########------- rrtf1 10m/0m or wt10m vs wt 0m -------------------------------------------------------------######
#Note that we subtract con because we defined it as reference and the change is wrt con
#we will use alpha=0.05 because that is the pvalue cutoff we will use later
#We will use IHW correction to improve power
library("IHW")
res_wt10vswt0 <- results(dds_keep, contrast = wt10 - wt0, parallel=TRUE, alpha=0.05,filterFun=ihw)
summary(res_wt10vswt0)
res_wt10vswt0$log2FoldChange[rownames(res_wt10vswt0) == "AT4G34410"]

res_rrtf0vswt0 <- results(dds_keep, contrast = rrtf0 - wt0, parallel=TRUE, alpha=0.05,filterFun=ihw)
summary(res_rrtf0vswt0)
res_rrtf0vswt0$log2FoldChange[rownames(res_rrtf0vswt0) == "AT4G34410"]

res_rrtf10vswt0 <- results(dds_keep, contrast = rrtf10 - wt0, parallel=TRUE, alpha=0.05,filterFun=ihw)
summary(res_rrtf10vswt0)
res_rrtf10vswt0$log2FoldChange[rownames(res_rrtf10vswt0) == "AT4G34410"]

res_rrtf10vsrrtf0 <- results(dds_keep, contrast = rrtf10 - rrtf0, parallel=TRUE, alpha=0.05,filterFun=ihw)
summary(res_rrtf10vsrrtf0)
res_rrtf10vsrrtf0$log2FoldChange[rownames(res_rrtf10vsrrtf0) == "AT4G34410"]

res_rrtf10vswt10 <- results(dds_keep, contrast = rrtf10 - wt10, parallel=TRUE, alpha=0.05,filterFun=ihw)
summary(res_rrtf10vswt10)
res_rrtf10vswt10$log2FoldChange[rownames(res_rrtf10vswt10) == "AT4G34410"]

#Filter genes with padj<0.05 and log2FC > 1
res_wt10vswt0 <- subset(res_wt10vswt0, padj < 0.05 & abs(log2FoldChange) > 1)
res_rrtf10vsrrtf0 <- subset(res_rrtf10vsrrtf0, padj < 0.05 & abs(log2FoldChange) > 1)
res_rrtf0vswt0 <- subset(res_rrtf0vswt0, padj < 0.05 & abs(log2FoldChange) > 1)
res_rrtf10vswt10 <- as.data.frame(subset(res_rrtf10vswt10, padj < 0.05 & abs(log2FoldChange) > 1))

write.csv(res_wt10vswt0, file = "res_wt10vswt0.csv")
write.csv(res_rrtf10vsrrtf0, file = "res_rrtf10vsrrtf0.csv")
write.csv(res_rrtf0vswt0, file = "res_rrtf0vswt0.csv")
res_rrtf0vswt0 <- as.data.frame(res_rrtf0vswt0)
res_wt10vswt0 <- as.data.frame(res_wt10vswt0)
res_rrtf10vsrrtf0 <- as.data.frame(res_rrtf10vsrrtf0)
res_rrtf10vsrrtf0 <- as.data.frame(res_rrtf10vsrrtf0)
~~~~~~~~~~~~~~~~~~~~~
  library("biomaRt")
ensembl <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org")
listAttributes(ensembl)
gene_info <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'description'), 
  filters = 'ensembl_gene_id', 
  values = rownames(counts), 
  mart = ensembl
)

res_wt10vswt0 <- res_wt10vswt0[,c(2,6)]
res_wt10vswt0$ensembl_gene_id <- ""
res_wt10vswt0$ensembl_gene_id <- rownames(res_wt10vswt0)
res_wt10vswt0_id <- merge(res_wt10vswt0, gene_info, by.x = "ensembl_gene_id", all.x = T)
write.csv(res_wt10vswt0_id, file = "res_wt10vswt0.csv")

res_rrtf0vswt0 <- res_rrtf0vswt0[,c(2,6)]
res_rrtf0vswt0$ensembl_gene_id <- ""
res_rrtf0vswt0$ensembl_gene_id <- rownames(res_rrtf0vswt0)
res_rrtf0vswt0_id <- merge(res_rrtf0vswt0, gene_info, by.x = "ensembl_gene_id", all.x = T)
write.csv(res_rrtf0vswt0_id, file = "res_rrtf0vswt0.csv")

res_rrtf10vsrrtf0 <- res_rrtf10vsrrtf0[,c(2,6)]
res_rrtf10vsrrtf0$ensembl_gene_id <- ""
res_rrtf10vsrrtf0$ensembl_gene_id <- rownames(res_rrtf10vsrrtf0)
res_rrtf10vsrrtf0_id <- merge(res_rrtf10vsrrtf0, gene_info, by.x = "ensembl_gene_id", all.x = T)
write.csv(res_rrtf10vsrrtf0_id, file = "res_rrtf10vsrrtf0.csv")

res_touch10 <- read.csv(file="res_touch10.csv", sep=",", header=TRUE,row.names = 1)[,c(2,6)]

res_rrtf10vswt10 <- res_rrtf10vswt10[,c(2,6)]
res_rrtf10vswt10$ensembl_gene_id <- ""
res_rrtf10vswt10$ensembl_gene_id <- rownames(res_rrtf10vswt10)
res_rrtf10vswt10_id <- merge(res_rrtf10vswt10, gene_info, by.x = "ensembl_gene_id", all.x = T)
write.csv(res_rrtf10vswt10_id, file = "res_rrtf10vswt10.csv")


#Obtain normalized counts from DESeq2
normalized_counts <- counts(dds_keep, normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts$ensembl_gene_id <- ""
normalized_counts$ensembl_gene_id <- rownames(normalized_counts)
normalized_counts <- merge(normalized_counts, gene_info, by.x = "ensembl_gene_id", all.x = T)
write.csv(normalized_counts, "normalized_counts.csv")

#####Merge publicaly available RNA-seq datasets####
library(dplyr)
data_frames <- list(res_touch5, res_touch10, res_touch20, res_touch25,res_touch30,res_touch40,res_touch60,res_touch180)
filtered_data_frames <- lapply(data_frames, function(df) {
  df %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1)
})
#double-check
t5f <- subset(res_touch5, padj < 0.05 & abs(log2FoldChange)>1)
t10f <- subset(res_touch10, padj < 0.05 & abs(log2FoldChange)>1)

# Step 2: Merge all filtered data frames by gene names (row names)
# First, make sure each data frame has the row names as a column for merging
filtered_data_frames <- lapply(filtered_data_frames, function(df) {
  df <- df %>% mutate(gene = rownames(df))
  return(df)
})
library(purrr)
# Step 3: Merge all filtered data frames on the 'gene' column
merged_df <- reduce(filtered_data_frames, full_join, by = "gene")
# Optionally, set "gene" as row names and remove it as a column
rownames(merged_df) <- merged_df$gene
merged_df$gene <- NULL
colnames(merged_df) <- c(
  "log2FC_5", "padj_5", "log2FC_10", "padj_10", "log2FC_20", "padj_20",
  "log2FC_25", "padj_25", "log2FC_30", "padj_30", "log2FC_40", "padj_40",
  "log2FC_60", "padj_60", "log2FC_180", "padj_180"
)
write.csv(merged_df, file="merged_df.csv")

#Obtain normalized counts from DESeq2
normalized_counts <- counts(dds_keep, normalized=TRUE)
#Reduce the gene size that matches merged_df
filtered_counts <- subset(counts_merge, rownames(counts_merge) %in% rownames(merged_df))
write.csv(counts_merge, file="raw_counts_all.csv")

#######Check batch effects of merged files
counts <- read.table("counts_xu_at_q30.txt",header=TRUE,row.names = 1)[,-c(1:5)]
setwd("~/../2023_Wang_RNAseq")
counts1 <- read.table("counts_combine_q30.txt",header=TRUE,row.names = 1)[,-c(1:5)]
setwd("~/../2019_Xu")
at_metadata <- read.table("at_metaData_merge.csv",header=TRUE,row.names = 1,sep=',')
at_metadata$time <- factor(at_metadata$time)
at_metadata$batch <- factor(at_metadata$batch)
# Check reference factor levels, if not 'con' for treatment and '1' for time relevel
print(at_metadata$time)
at_metadata$time <- relevel(at_metadata$time, ref = '0')
# next three steps should be in exact order 
sample_order <- order(at_metadata$time)
at_metadata <- at_metadata[sample_order, , drop=F] #drop=F ensures that even if you're subsetting a single row, the result will still be a data frame,not a list.
#merge two counts' files
counts_merge <- merge(counts, counts1, by="row.names", all=T)
rownames(counts_merge) <- counts_merge[,1]
counts_merge <- counts_merge[,-1]
counts_merge <- counts_merge[, rownames(at_metadata)]
colnames(counts_merge)
rownames(at_metadata)
all(colnames(counts_merge) == rownames(at_metadata))

dds <- DESeqDataSetFromMatrix(countData = counts_merge,
                              colData = at_metadata,
                              design = ~ time + batch)
keep <- rowSums(counts(dds) > 10) >= 3
dds_keep <- dds[keep,]
#Normalized counts
dds_keep <- DESeq(dds_keep,betaPrior=FALSE, parallel=TRUE)
resultsNames(dds_keep)
# Check sample outlier using Cooks distance
summary(results(dds_keep, alpha=0.05))
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_keep)[["cooks"]]), range=0, las=2)

#remove batcheffect
library(limma)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
y<-plotPCA(vsd, intgroup = c("batch","time"))
y + geom_text(aes(label=substr(name, start = 1, stop = 12)),vjust=2,check_overlap = T,size = 2)

####TPM calculation
counts1 <- read.table("counts_rrtf_at_q30.txt",header=TRUE,row.names = 1)
x <- counts.mat / gene.length

#####significant genes from DESeq2 all constrasts, padj <0.05####
sig_w10_w0 <- as.data.frame(subset(res_wt10vswt0, padj < 0.05))
sig_r10_r0 <- as.data.frame(subset(res_rrtf10vsrrtf0, padj < 0.05))
sig_r0_w0 <- as.data.frame(subset(res_rrtf0vswt0, padj < 0.05))
sig_r10_w10 <- as.data.frame(subset(res_rrtf10vswt10, padj < 0.05))

sig_w10_w0 <- sig_w10_w0[,c(2,6)]
sig_w10_w0$ensembl_gene_id <- ""
sig_w10_w0$ensembl_gene_id <- rownames(sig_w10_w0)
sig_r10_r0 <- sig_r10_r0[,c(2,6)]
sig_r10_r0$ensembl_gene_id <- ""
sig_r10_r0$ensembl_gene_id <- rownames(sig_r10_r0)
sig_r0_w0 <- sig_r0_w0[,c(2,6)]
sig_r0_w0$ensembl_gene_id <- ""
sig_r0_w0$ensembl_gene_id <- rownames(sig_r0_w0)
sig_r10_w10 <- sig_r10_w10[,c(2,6)]
sig_r10_w10$ensembl_gene_id <- ""
sig_r10_w10$ensembl_gene_id <- rownames(sig_r10_w10)


sig_1 <- merge(sig_w10_w0, sig_r10_r0, by = "row.names", all.x = T, all.y = T)
rownames(sig_1) <- sig_1[,1]
colnames(sig_1) <- c("f","log2fc_w10_w0","padj_w10_w0","fff","log2fc_r10_r0","padj_r10_r0","fffff")
sig_1 <- sig_1[,c(2,3,5,6)]
sig_2 <- merge(sig_1, sig_r0_w0, by = "row.names", all.x = T, all.y = T)
rownames(sig_2) <- sig_2[,1]
sig_2 <- sig_2[,2:7]
colnames(sig_2)[c(5,6)] <- c("log2fc_r0_w0","padj_r0_w0")
sig_3 <- merge(sig_2, sig_r10_w10, by = "row.names", all.x = T, all.y = T)
rownames(sig_3) <- sig_3[,1]
sig_3 <- sig_3[,2:9]
colnames(sig_3)[c(7,8)] <- c("log2fc_r10_w10","padj_r10_w10")
sig_3 <- sig_3[,-1]

sig_3$ensembl_gene_id <- ""
sig_3$ensembl_gene_id <- rownames(sig_3)
sig_3 <- merge(sig_3, gene_info, by.x = "ensembl_gene_id", all.x = T)

keep_sig <- rownames(counts) %in% rownames(sig_3)
dds_keep_sig <- dds[keep_sig,]
#Normalized counts
dds_keep_sig <- DESeq(dds_keep_sig,betaPrior=FALSE, parallel=TRUE)
resultsNames(dds_keep_sig)
normalized_counts_sig <- counts(dds_keep_sig, normalized=TRUE)
normalized_counts_sig <- as.data.frame(normalized_counts_sig)
