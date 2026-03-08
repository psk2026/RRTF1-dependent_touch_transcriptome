library(edgeR)
###First, do check raw count files.####
counts <- read.table("counts_rrtf_at_q30.txt",header=TRUE,row.names = 1)[,-c(1:5)]
at_metadata <- read.table("at_metaData.csv",header=TRUE,row.names = 1,sep=',')
at_metadata$merge <- factor(at_metadata$merge)

keep <- rowSums(counts(dds) > 10) >= 3

dge <- DGEList(counts = counts, group = factor(at_metadata$merge))
#Filter out the genes with low counts
keep <- filterByExpr(y = dge, group=factor(at_metadata$merge))
dge <- dge[keep, , keep.lib.sizes=FALSE]

#Normalization and effective library sizesP
dge <- calcNormFactors(object = dge)
dge$samples
plotMDS(dge)

design <- model.matrix(~ 0 + factor(at_metadata$merge))
colnames(design) <- levels(factor(at_metadata$merge))

#Model fitting and estimating dispersions
dge <- estimateDisp(y = dge, design, robust = T)
dge$common.dispersion
plotBCV(dge)

#edgeR-glm pipeline
fit <- glmQLFit(dge, design, robust=TRUE)
fit$dispersion
head(fit$coefficients)
plotQLDisp(fit)

con <- makeContrasts(wt10 - wt0, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))
qlf$table[rownames(qlf$table) == "AT2G37430",]

tr <- glmTreat(fit, contrast=con, lfc=log2(1.3))
summary(decideTests(tr))
tr$table[rownames(tr$table) == "AT2G37430",]
plotMD(tr)
edge_wt10_wt0 <- topTags(tr, n=Inf)$table
edge_wt10_wt0 <- edge_wt10_wt0[edge_wt10_wt0$FDR < 0.05,]
edge_wt10_wt0 <- edge_wt10_wt0[abs(edge_wt10_wt0$unshrunk.logFC)>1,]
write.csv(edge_wt10_wt0, file="edge_wt10_wt0.csv")


con <- makeContrasts(rrtf10 - rrtf0, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))
qlf$table[rownames(qlf$table) == "AT2G37430",]
tr <- glmTreat(fit, contrast=con, lfc=log2(1.3))
summary(decideTests(tr))
tr$table[rownames(tr$table) == "AT2G37430",]
plotMD(tr)
edge_rrtf10_rrtf0 <- topTags(tr, n=Inf)$table
edge_rrtf10_rrtf0 <- edge_rrtf10_rrtf0[edge_rrtf10_rrtf0$FDR < 0.05,]
edge_rrtf10_rrtf0 <- edge_rrtf10_rrtf0[edge_rrtf10_rrtf0$FDR < 0.05,]
edge_rrtf10_rrtf0 <- edge_rrtf10_rrtf0[abs(edge_rrtf10_rrtf0$unshrunk.logFC)>1,]

con <- makeContrasts(rrtf10 - wt0, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))
qlf$table[rownames(qlf$table) == "AT2G37430",]
tr <- glmTreat(fit, contrast=con, lfc=log2(1.3))
summary(decideTests(tr))
tr$table[rownames(tr$table) == "AT2G37430",]
plotMD(tr)
edge_rrtf10_wt0 <- topTags(tr, n=Inf)$table
edge_rrtf10_wt0 <- edge_rrtf10_wt0[edge_rrtf10_wt0$FDR < 0.05,]
edge_rrtf10_wt0 <- edge_rrtf10_wt0[abs(edge_rrtf10_wt0$unshrunk.logFC)>1,]
write.csv(edge_rrtf10_wt0, file="edge_rrtf10_wt0.csv")

con <- makeContrasts(rrtf10 - wt10, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))
qlf$table[rownames(qlf$table) == "AT2G37430",]

tr <- glmTreat(fit, contrast=con, lfc=log2(1))
summary(decideTests(tr))
tr$table[rownames(tr$table) == "AT2G37430",]
plotMD(tr)
edge_rrtf10_wt10 <- topTags(tr, n=Inf)$table
edge_rrtf10_wt10 <- edge_rrtf10_wt10[edge_rrtf10_wt10$FDR < 0.05,]
edge_rrtf10_wt10 <- edge_rrtf10_wt10[abs(edge_rrtf10_wt10$logFC)>1,]
write.csv(edge_rrtf10_wt10, file="edge_rrtf10_wt10.csv")

con <- makeContrasts(rrtf0 - wt0, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))
qlf$table[rownames(qlf$table) == "AT2G37430",]

tr <- glmTreat(fit, contrast=con, lfc=log2(1.3))
summary(decideTests(tr))
tr$table[rownames(tr$table) == "AT2G37430",]
plotMD(tr)
edge_rrtf0_wt0 <- topTags(tr, n=Inf)$table
edge_rrtf0_wt0 <- edge_rrtf0_wt0[edge_rrtf0_wt0$FDR < 0.05,]
edge_rrtf0_wt0 <- edge_rrtf0_wt0[abs(edge_rrtf0_wt0$unshrunk.logFC)>1,]
write.csv(edge_rrtf0_wt0, file="edge_rrtf0_wt0.csv")

edge_wt10vswt0_rrtf10vswt0 <- merge(edge_wt10_wt0, edge_rrtf10_wt0, by="row.names", all=T)
colnames(edge_wt10vswt0_rrtf10vswt0) <- c("ensembl_gene_id", "log2FC_wt10_wt0", "unshrunk_log2FC_wt10_wt0", "logCPM_wt10_wt0", "pval_wt10_wt0", "FDR_wt10_wt0", "log2FC_rrtf10_wt0", "unshrunk_log2FC_rrtf10_wt0","logCPM_rrtf10_wt0", "pval_rrtf10_wt0", "FDR_rrtf10_wt0")
write.csv(edge_wt10vswt0_rrtf10vswt0, file="edge_wt10vswt0_rrtf10vswt0.csv")

merge_edge_geno <- merge(edge_wt10_wt0, edge_rrtf10_rrtf0, by= "row.names", all.x=T, all.y = T)
colnames(merge_edge_geno) <- c("ensembl_gene_id", "log2FC_wt10_wt0", "unshrunk_log2FC_wt10_wt0", "logCPM_wt10_wt0", "pval_wt10_wt0", "FDR_wt10_wt0", "log2FC_rrtf10_rrtf0", "unshrunk_log2FC_rrtf10_rrtf0","logCPM_rrtf10_rrtf0", "pval_rrtf10_rrtf0", "FDR_rrtf10_rrtf0")

merge_edge_geno <- merge(merge_edge_geno, gene_info, by.x = "ensembl_gene_id", all.x = T)
write.csv(merge_edge_geno1, "merge_edge_geno.csv")


#ANOVA-like testing
con <- makeContrasts(wt10_wt0 = wt10 - wt0, rrtf10_rrtf0 = rrtf10 - rrtf0, levels=design)
anov <- glmQLFTest(fit, contrast=con)
topTags(anov)
anov$table[rownames(anov$table) == "AT2G37430",]
write.csv(anov, "anov.csv")

#core touch-responsive genes > extract DEG from each DEG analysis
rrtf10vswt10_rrtfonly <- read.csv(file="rrtf10vswt10_rrtfonly.csv", sep=",")
core_VM <- read.csv(file="core_touch_genes_VM.csv",sep=",")

frrtf10 <- edge_rrtf10_wt0
frrtf10$gene <- ""
frrtf10$gene <- rownames(frrtf10)
frrtf10 <- frrtf10[frrtf10$FDR < 0.05,]

fwt10 <- edge_wt10_wt0
fwt10$gene <- ""
fwt10$gene <- rownames(fwt10)
fwt10 <- fwt10[fwt10$FDR < 0.05,]

frrtf0 <- edge_rrtf0_wt0
frrtf0$gene <- ""
frrtf0$gene <- rownames(frrtf0)
frrtf0 <- frrtf0[frrtf0$FDR < 0.05,]

rrtf10vswt10_rrtfonly_rrtf10 <- merge(frrtf10, core_VM, by="gene") #overlap with rrtf10(vswt0) and core_touch_VM
rrtf10vswt10_rrtfonly_rrtf10 <- rrtf10vswt10_rrtfonly_rrtf10[c(1,2)]
colnames(rrtf10vswt10_rrtfonly_rrtf10)[2] <- c("rrtf10")

rrtf10vswt10_rrtfonly_wt10 <- merge(fwt10, core_VM, by = "gene")
rrtf10vswt10_rrtfonly_wt10 <- rrtf10vswt10_rrtfonly_wt10[c(1,2)]
colnames(rrtf10vswt10_rrtfonly_wt10)[2] <- c("wt10")

rrtf10vswt10_rrtfonly_rrtf0 <- merge(frrtf0, core_VM, by = "gene")
rrtf10vswt10_rrtfonly_rrtf0 <- rrtf10vswt10_rrtfonly_rrtf0[c(1,2)]
colnames(rrtf10vswt10_rrtfonly_rrtf0)[2] <- c("rrtf0")

write.csv(rrtf10vswt10_rrtfonly_rrtf10, "rrtf10vswt10_rrtfonly_rrtf10.csv")
write.csv(rrtf10vswt10_rrtfonly_wt10, "rrtf10vswt10_rrtfonly_wt10.csv")
write.csv(rrtf10vswt10_rrtfonly_rrtf0, "rrtf10vswt10_rrtfonly_rrtf0.csv")

core_wt10_rrtf10_rrtf0 <- merge(rrtf10vswt10_rrtfonly_wt10, rrtf10vswt10_rrtfonly_rrtf10, by = "gene", all=T)
core_wt10_rrtf10_rrtf0 <- merge(core_wt10_rrtf10_rrtf0, rrtf10vswt10_rrtfonly_rrtf0, by = "gene", all=T)
write.csv(core_wt10_rrtf10_rrtf0, "core_wt10_rrtf10_rrtf0.csv")

###non core touch-responsive genes
filtered_genes_vector <- all_my_genes_vector[!all_my_genes_vector %in% genes_to_exclude_list]
frrtf10 <- edge_rrtf10_wt0
frrtf10$gene <- ""
frrtf10$gene <- rownames(frrtf10)
frrtf10 <- frrtf10[frrtf10$FDR < 0.05,]

fwt10 <- edge_wt10_wt0
fwt10$gene <- ""
fwt10$gene <- rownames(fwt10)
fwt10 <- fwt10[fwt10$FDR < 0.05,]

frrtf0 <- edge_rrtf0_wt0
frrtf0$gene <- ""
frrtf0$gene <- rownames(frrtf0)
frrtf0 <- frrtf0[frrtf0$FDR < 0.05,]

library(dplyr)
#overlap with rrtf10(vswt0) and core_touch_VM
rrtf10vswt10_rrtfonly_rrtf10 <- dplyr::anti_join(frrtf10, core_VM, by = "gene")
rrtf10vswt10_rrtfonly_rrtf10 <- rrtf10vswt10_rrtfonly_rrtf10[c(1,6)]
colnames(rrtf10vswt10_rrtfonly_rrtf10)[1] <- c("rrtf10")

rrtf10vswt10_rrtfonly_wt10 <- dplyr::anti_join(fwt10, core_VM, by = "gene")
rrtf10vswt10_rrtfonly_wt10 <- rrtf10vswt10_rrtfonly_wt10[c(1,6)]
colnames(rrtf10vswt10_rrtfonly_wt10)[1] <- c("wt10")

rrtf10vswt10_rrtfonly_rrtf0 <- dplyr::anti_join(frrtf0, core_VM, by = "gene")
rrtf10vswt10_rrtfonly_rrtf0 <- rrtf10vswt10_rrtfonly_rrtf0[c(1,6)]
colnames(rrtf10vswt10_rrtfonly_rrtf0)[1] <- c("rrtf0")

write.csv(rrtf10vswt10_rrtfonly_rrtf10, "rrtf10vswt10_rrtfonly_rrtf10.csv")
write.csv(rrtf10vswt10_rrtfonly_wt10, "rrtf10vswt10_rrtfonly_wt10.csv")
write.csv(rrtf10vswt10_rrtfonly_rrtf0, "rrtf10vswt10_rrtfonly_rrtf0.csv")

noncore_wt10_rrtf10_rrtf0 <- merge(rrtf10vswt10_rrtfonly_wt10, rrtf10vswt10_rrtfonly_rrtf10, by = "gene", all=T)
noncore_wt10_rrtf10_rrtf0 <- merge(noncore_wt10_rrtf10_rrtf0, rrtf10vswt10_rrtfonly_rrtf0, by = "gene", all=T)
write.csv(noncore_wt10_rrtf10_rrtf0, "noncore_wt10_rrtf10_rrtf0.csv")


edgeDEGmerge_vswt0 <- merge(edge_wt10_wt0, edge_rrtf10_wt0, by="row.names", all=T)
edgeDEGmerge_vswt0 <- edgeDEGmerge_vswt0[,c(1,2,3,6,7,8,11)]
rownames(edgeDEGmerge_vswt0) <- edgeDEGmerge_vswt0[,1]
edgeDEGmerge_vswt0 <- edgeDEGmerge_vswt0[,-1]
colnames(edgeDEGmerge_vswt0) <- c("logFC_wt10", "unshrunk.logFC_wt10", "FDR_wt10", "logFC_rrtf10","unshrunk.logFC_rrtf10","FDR_rrtf10")
edgeDEGmerge_vswt0 <- merge(edgeDEGmerge_vswt0, edge_rrtf0_wt0, by="row.names", all=T)
write.csv(edgeDEGmerge_vswt0, "edgeDEGmerge_vswt0.csv")
colnames(edgeDEGmerge_vswt0)
