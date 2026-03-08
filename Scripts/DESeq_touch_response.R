#####DESeq_for_Merged_Bam_file
setwd("~/../Thigmo_exp/2025_RNA_seq/Merge_bam/")

#columns 2:6 from data_counts contain values other than counts
at_counts <- read.table("counts_merged_at_q30.txt",header=TRUE,row.names = 1)[,-c(1:5)]
metadata <- read.table("at_metaData_onefactor.csv",header=TRUE,row.names = 1,sep=',')
metadata$time <- factor(metadata$time)

# Check reference factor levels, if not 'con' for treatment
print(metadata$time)
metadata$time <- relevel(metadata$time, ref = "con")

metadata <- droplevels(metadata)
print(metadata)
# colnames should be in exact order 
at_counts <- at_counts[,rownames(metadata)]

library(dplyr)
library(ggplot2)
library(DESeq2)
# Use multiple cores for DEseq2
library("BiocParallel")
register(MulticoreParam(8))

dds <- DESeqDataSetFromMatrix(countData = at_counts,
                              colData = metadata,
                              design = ~ time)

# At least 10 counts in 3 samples
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
plotPCA(vsd, intgroup=c("time"))
y<-plotPCA(vsd, intgroup=c("time"))
y + geom_text(aes(label=substr(name, start = 1, stop = 6)),vjust=2,check_overlap = T,size = 2)
cbind(y$data$time,y$data$PC1,y$data$PC2)

### Using all samples / get the model matrix
mod_mat <- model.matrix(design(dds_keep), colData(dds_keep))

# Define coefficient vectors for each condition
con <- colMeans(mod_mat[dds_keep$time == "con", ])
t5min <- colMeans(mod_mat[dds_keep$time == "t1_5min", ])
t10min <- colMeans(mod_mat[dds_keep$time == "t2_10min", ])
t30min <- colMeans(mod_mat[dds_keep$time == "t3_30min", ])
########------- touch vs con -------------------------------------------------------------######
#Note that we subtract con because we defined it as reference.
#we will use alpha=0.05 because that is the pvalue cutoff we will use later.
#We will use IHW correction to improve power.
library("IHW")
res_5min <- results(dds_keep, contrast = t5min - con, parallel=TRUE, alpha=0.05,filterFun=ihw)
summary(res_5min)
res_10min <- results(dds_keep, contrast = t10min - con, parallel=TRUE, alpha=0.05,filterFun=ihw)
summary(res_10min)
res_30min <- results(dds_keep, contrast = t30min - con, parallel=TRUE, alpha=0.05,filterFun=ihw)
summary(res_30min)

#Number of genes below padj<0.05
sum(res_5min$padj < 0.05, na.rm=TRUE)
sum(res_10min$padj < 0.05, na.rm=TRUE)
sum(res_30min$padj < 0.05, na.rm=TRUE)

#order the results using adj pvalues and log2FC.
res_5min_filter <- subset(res_5min[,c(2,3,6)], padj < 0.05 & abs(log2FoldChange)>1)
res_5min_order <- res_5min_filter[order(res_5min_filter$log2FoldChange),]

res_10min_filter <- subset(res_10min[,c(2,3,6)], padj < 0.05 & abs(log2FoldChange)>1)
res_10min_order <- res_10min_filter[order(res_10min_filter$log2FoldChange),]

res_30min_filter <- subset(res_30min[,c(2,3,6)], padj < 0.05 & abs(log2FoldChange)>1)
res_30min_order <- res_30min_filter[order(res_30min_filter$log2FoldChange),]

res_5min_order <- as.data.frame(res_5min_order)
res_10min_order <- as.data.frame(res_10min_order)
res_30min_order <- as.data.frame(res_30min_order)

sum(res_5min_order$log2FoldChange > 0)
sum(res_5min_order$log2FoldChange < 0)

sum(res_10min_order$log2FoldChange > 0)
sum(res_10min_order$log2FoldChange < 0)

sum(res_30min_order$log2FoldChange > 0)
sum(res_30min_order$log2FoldChange < 0)
write.csv(res_5min_order, file="res_5min_order.csv")
write.csv(res_10min_order, file="res_10min_order.csv")
write.csv(res_30min_order, file="res_30min_order.csv")

core_touch_genes_VM <- read.table("core_touch_genes_VM.csv",header=TRUE, sep=",")
rownames(core_touch_genes_VM) <- core_touch_genes_VM[,1]
overlap5min <- merge(res_5min_order, core_touch_genes_VM, by = "row.names") #195/218
overlap10min <- merge(core_touch_genes_VM, res_10min_order, by = "row.names") #541/673
overlap30min <- merge(core_touch_genes_VM, res_30min_order, by = "row.names") #869/1417

#TF comparison
TFDB <- read.table("Ath_TF_list_TFDB.txt",header=T)
rownames(TFDB) <- TFDB[,2]
res_5min_TF <- subset(res_5min_order, rownames(res_5min_order) %in% TFDB$Gene_ID) #47/218 = 21.6%
res_10min_TF <- subset(res_10min_order, rownames(res_10min_order) %in% TFDB$Gene_ID) #96/673 = 14.2%
res_30min_TF <- subset(res_30min_order, rownames(res_30min_order) %in% TFDB$Gene_ID) #166/1417 = 11.7%

#Total TF fraction
total_deg <- merge(res_10min_order, res_30min_order, by = "row.names", all.x=T, all.y=T)
rownames(total_deg) <- total_deg[,1]
total_deg <- total_deg[,-1]
total_deg <- merge(res_5min_order, total_deg, by = "row.names", all.x=T, all.y=T)
rownames(total_deg) <- total_deg[,1]
total_deg <- total_deg[,-1]
total_TF <- subset(total_deg, rownames(total_deg) %in% TFDB$Gene_ID) #183/1574 = 11.6%

###normalized counts
norm_dds_keep <- estimateSizeFactors(dds_keep)
sizeFactors(norm_dds_keep)
normalized_counts <- counts(norm_dds_keep, normalized=T)
normalized_counts <- as.data.frame(normalized_counts)
deg_counts <- normalized_counts %>% filter(rownames(normalized_counts) %in% rownames(total_deg))
