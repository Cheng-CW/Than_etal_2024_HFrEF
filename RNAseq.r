# Than April 2024 Male versus Female

setwd("~/Desktop/")
library(DESeq2)
library(dplyr)
countData = read.csv("readcounts.txt", header=TRUE, sep="\t", row.names = 1, check.names=FALSE)
sampleCondition = read.csv("metaData_RNAseq_male_vs_female.csv", header=TRUE)
dim(countData)
dim(sampleCondition)
ids <- read.csv("id_name.txt", header=TRUE, sep="\t", row.names=1)

# Construct principal component analysis using raw counts across all samples.
library(ggfortify)
library(ggrepel)
pcDat <- prcomp(t(countData))
head(pcDat)
autoplot(pcDat, size=2, label = TRUE, label.repel=T) 

# Identify DEGs
sample_analysis_1_DE <- DESeqDataSetFromMatrix(countData=countData, colData=sampleCondition, design= ~Condition)
as.data.frame(colData(sample_analysis_1_DE))
analysis_1 <- sample_analysis_1_DE[rowSums(DESeq2::counts(sample_analysis_1_DE)) >= 10, ]
analysis_1$Condition = relevel(analysis_1$Condition, ref = "female_CON")
analysis_1_dds <- DESeq(analysis_1)
resultsNames(analysis_1_dds)
colSums(assay(analysis_1_dds))
colData(analysis_1_dds)
normalised_analysis_1 <- counts(analysis_1_dds, normalized=TRUE) 
head(normalised_analysis_1)
analysis_1_DEG <-results(analysis_1_dds)
analysis_1_DEGs <- analysis_1_DEG[order(analysis_1_DEG$padj),]
analysis_1_result <- as.data.frame(dplyr::mutate(as.data.frame(analysis_1_DEGs ), sig=ifelse(analysis_1_DEGs$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(analysis_1_DEGs))
sum(analysis_1_result$padj<0.05, na.rm = T)
sum(analysis_1_result$padj<0.05 & analysis_1_result$log2FoldChange > 0, na.rm = T)
sum(analysis_1_result$padj<0.05 & analysis_1_result$log2FoldChange < 0, na.rm = T)
summary(analysis_1_result)

###### Shrinkage of log2fold change
resLFC  <- lfcShrink(analysis_1_dds, coef= "Condition_M_vs_F", type="apeglm")
head(resLFC)
plotMA(resLFC)
sum(resLFC$padj<0.05, na.rm = T)
resLFC_DEGs <- resLFC[order(resLFC$padj),]
head(resLFC_DEGs)

# PCA
analysis_1_PCA <- vst(analysis_1_dds, blind=FALSE)

pca_1 <- plotPCA(analysis_1_PCA, intgroup=c("Condition"), returnData=TRUE) 

percentVar_1 <- round(100 * attr(pca_1, "percentVar"))

# Batch and name
ggplot(pca_1, aes(PC1, PC2, color=Condition)) + 
  geom_text_repel(aes(label=name)) + 
  theme(text = element_text(size = 10)) +
  geom_point(size = 5) + ggtitle("PCA") +
  xlab(paste0("PC1: ",percentVar_1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_1[2],"% variance"))

# Batch 
ggplot(pca_1, aes(PC1, PC2, color=Condition)) + 
  theme(text = element_text(size = 5)) +
  geom_point(size = 5) + ggtitle("PCA") +
  xlab(paste0("PC1: ",percentVar_1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_1[2],"% variance")) + 
  coord_fixed(ratio=1)

# Condition
ggplot(pca_1, aes(PC1, PC2, color=Condition)) + 
  theme(text = element_text(size = 5)) +
  geom_point(size = 5) + ggtitle("PCA") +
  xlab(paste0("PC1: ",percentVar_1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_1[2],"% variance")) + 
  coord_fixed(ratio=2)

# PLS-DA
library(mixOmics)
normalised_analysis_1_transpose <- t(normalised_analysis_1)
plsda_analysis_1 = normalised_analysis_1_transpose
outcome_analysis_1 <- sampleCondition$Condition
plsda_res_analysis_1 = plsda(plsda_analysis_1, outcome_analysis_1, ncomp = 3, max.iter = 1000)
plotIndiv(plsda_res_analysis_1, ind.names = TRUE, legend = TRUE, ellipse = TRUE, title = 'PLS-DA')

# Enrichment analysis
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(ensembldb)
library(ReactomePA)
DEGs_analysis_1 <- subset(analysis_1_result, padj < 0.05)
dim(DEGs_analysis_1)
ensembl_analysis_1 <- rownames(DEGs_analysis_1)
ensembl_analysis_1 <- sub("*\\..*", "", ensembl_analysis_1)
ensembl_analysis_1
id_analysis_1 <- bitr(ensembl_analysis_1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
Enrichment_DEGs_analysis_1 <- enrichGO(gene = as.vector(id_analysis_1$ENTREZID), 
                                       keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", 
                                       pvalueCutoff=0.05, qvalueCutoff=0.05, minGSSize=2, readable = TRUE)
cluster_summary_analysis_1 <- data.frame(Enrichment_DEGs_analysis_1)
dim(cluster_summary_analysis_1)
ego_analysis_1 <- setReadable(Enrichment_DEGs_analysis_1, OrgDb = org.Hs.eg.db)
head(ego_analysis_1)
clusterProfiler::dotplot(ego_analysis_1, showCategory=20)
cnetplot(ego_analysis_1)

clusterProfiler::dotplot(ego_analysis_1, showCategory=10, font.size=23) +
  scale_y_discrete(labels = function(y) lapply(strwrap(y, width = 30, simplify = FALSE),  paste, collapse="\n")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

write.csv(ego_analysis_1, file="enrichment.csv")

# KEGG
kegg_organism <- "hsa"
kegg_analysis_1 <- bitr(ensembl_analysis_1, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
kegg_analysis_1 <- kegg_analysis_1[!duplicated(kegg_analysis_1[c("ENSEMBL")]),]
DEGs_analysis_1$ensembl <- rownames(DEGs_analysis_1)
DEGs_analysis_1$ensembl <- sub("*\\..*", "", DEGs_analysis_1$ensembl)
row.names(DEGs_analysis_1) <- DEGs_analysis_1$ensembl
DEGs_analysis_1$ensembl <- NULL
kegg_analysis_1_a <- DEGs_analysis_1[rownames(DEGs_analysis_1) %in% kegg_analysis_1$ENSEMBL,]
kegg_analysis_1_a$Y <- kegg_analysis_1$ENTREZID
FC_analysis_1 <- kegg_analysis_1_a$log2FoldChange
names(FC_analysis_1) <- kegg_analysis_1_a$Y
kegg_gene_list_analysis_1 <- sort(FC_analysis_1, decreasing = TRUE)
kegg_analysis_1 <- gseKEGG(geneList = kegg_gene_list_analysis_1,
                           organism     = kegg_organism,
                           nPerm        = 1000,
                           minGSSize    = 2,
                           maxGSSize    = 800,
                           pvalueCutoff = 0.10,
                           pAdjustMethod = "none",
                           keyType       = "ncbi-geneid")
gseaKEGG_results_analysis_1 <- kegg_analysis_1@result
head(gseaKEGG_results_analysis_1)
kegg_readable_analysis_1 <- setReadable(kegg_analysis_1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_cluster_summary_analysis_1 <- data.frame(kegg_readable_analysis_1)
clusterProfiler::dotplot(kegg_readable_analysis_1, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
write.csv(kegg_cluster_summary_analysis_1, file="KEGG.csv")

# REACTOME
reactome_analysis_1 <- enrichPathway(gene=kegg_analysis_1_a$Y, pvalueCutoff=0.05, pAdjustMethod="none", readable=T, organism="human")
reactome_analaysis_1_A <- (as.data.frame(reactome_analysis_1))
head(reactome_analysis_1)
dim(reactome_analysis_1)
clusterProfiler::dotplot(reactome_analysis_1, showCategory=10, font.size=23) +
  scale_y_discrete(labels = function(y) lapply(strwrap(y, width = 40, simplify = FALSE),  paste, collapse="\n")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
write.csv(reactome_analysis_1, file="reactome.csv")

# Organise data
analysis_1_merge <- merge(analysis_1_result, normalised_analysis_1, by=0, all=TRUE)
rownames(analysis_1_merge) = analysis_1_merge$Row.names
analysis_1_merge$Row.names <- NULL
head(analysis_1_merge)
analysis_1_merge_1 <- merge(analysis_1_merge, ids, by=0, all=TRUE)
rownames(analysis_1_merge_1) <- analysis_1_merge_1$Row.names
analysis_1_merge_1$Row.names <- NULL
analysis_1_merge_2 <- analysis_1_merge_1 %>% arrange(padj)
analysis_1_merge_3 <- na.omit(analysis_1_merge_2)
final_1 <- as.data.frame(dplyr::mutate(as.data.frame(analysis_1_merge_3), 
                                       M =ifelse(analysis_1_merge_3$log2FoldChange>0, "UP", "DOWN")), 
                         row.names=rownames(analysis_1_merge_3))
write.csv(final_1, file="DEGS.csv")

# Volcano plot 
volcano_plot_1 <- analysis_1_merge_3[ which(analysis_1_merge_3$padj < 0.05),]
rownames(volcano_plot_1) = volcano_plot_1$gene_symbol

# Duplicate names  rownames(volcano_plot_1) <- make.names(volcano_plot_1$gene_symbol, unique = TRUE) 

p <- ggplot(analysis_1_merge_3, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col = sig)) + 
  scale_color_manual(values = c("red", "black")) + 
  geom_vline(xintercept=c(-0,0), linetype="dotted") + 
  theme(legend.position = "none", text = element_text(size=20))
p + geom_text_repel(data=volcano_plot_1[1:30, ], aes(label=rownames(volcano_plot_1[1:30, ])))

#####################################################################################################################
############################################ COMBINING TWO DATASET #################################################
#####################################################################################################################

setwd("~/Desktop/")
library(DESeq2)
library(dplyr)
ids <- read.csv("id_name.txt", header=TRUE, sep="\t", row.names=1)
countData_1 = read.csv("counts_malevsfemale.csv", header=TRUE, row.names = 1, check.names=FALSE)
sampleCondition_1 = read.csv("sample_info_malevsfemale.csv", header=TRUE)
dim(countData_1)
dim(sampleCondition_1)

countData_2 = read.csv("readcounts.txt", header=TRUE, sep="\t", row.names = 1, check.names=FALSE)
sampleCondition_2 = read.csv("metaData_RNAseq_male_vs_female.csv", header=TRUE)
dim(countData_2)
dim(sampleCondition_2)

# Merge both datasets (count)
overall_counts <- merge(countData_1, countData_2, by=0, all=TRUE)
row.names(overall_counts) <- overall_counts$Row.names
overall_counts$Row.names <- NULL
dim(overall_counts)

# Merge both datasets (metaData)
overall_sample_conditions <- rbind(sampleCondition_1, sampleCondition_2)
overall_sample_conditions$Batch <- factor(overall_sample_conditions$Batch)
dim(overall_sample_conditions)
head(overall_sample_conditions)

# Construct principal component analysis using raw counts across all samples.
library(ggfortify)
library(ggrepel)
pcDat <- prcomp(t(overall_counts))
head(pcDat)
autoplot(pcDat, size=2, label = TRUE, label.repel=T) 

# Remove batch effect using ComBat-seq
library(sva)

##############################
# Analysis 1 (Male_HF VERSUS Male_CON) 
##############################
sample_analysis_1 <- overall_sample_conditions[which(overall_sample_conditions$Condition == "male_HF" | overall_sample_conditions$Condition == "male_CON"),]
head(sample_analysis_1)
dim(sample_analysis_1)
table(sample_analysis_1$Condition)
sample_analysis_1_ids <- sample_analysis_1$Sample_File
sample_analysis_1_ids <- as.vector(t(sample_analysis_1_ids))
head(sample_analysis_1_ids)
sample_analysis_1_count <- overall_counts[, sample_analysis_1_ids]
head(sample_analysis_1_count)
dim(sample_analysis_1_count)

### Correct batch effect using ComBat-seq
#sample_analysis_1_count <- ComBat_seq(counts = sample_analysis_1_count, batch = sample_analysis_1$Batch)

# Identify DEGs
sample_analysis_1$Batch <- factor(sample_analysis_1$Batch)
sample_analysis_1_DE <- DESeqDataSetFromMatrix(countData=sample_analysis_1_count, colData=sample_analysis_1, design=~Condition)
as.data.frame(colData(sample_analysis_1_DE))
analysis_1 <- sample_analysis_1_DE[rowSums(DESeq2::counts(sample_analysis_1_DE)) > 10, ]
analysis_1_dds <- DESeq(analysis_1)
resultsNames(analysis_1_dds)
colSums(assay(analysis_1_dds))
colData(analysis_1_dds)
normalised_analysis_1 <- counts(analysis_1_dds, normalized=TRUE) 
head(normalised_analysis_1)
analysis_1_DEG <-results(analysis_1_dds)
analysis_1_DEGs <- analysis_1_DEG[order(analysis_1_DEG$padj),]
analysis_1_result <- as.data.frame(dplyr::mutate(as.data.frame(analysis_1_DEGs ), sig=ifelse(analysis_1_DEGs$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(analysis_1_DEGs))
sum(analysis_1_result$padj<0.05, na.rm = T)
summary(analysis_1_result)
sum(analysis_1_result$padj<0.05, na.rm = T)
sum(analysis_1_result$padj<0.05 & analysis_1_result$log2FoldChange > 0, na.rm = T)
sum(analysis_1_result$padj<0.05 & analysis_1_result$log2FoldChange < 0, na.rm = T)

# Correct batch effect using sva
dat  <- counts(analysis_1_dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ Condition, colData(analysis_1_dds))
mod0 <- model.matrix(~   1, colData(analysis_1_dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
svseq$sv
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ analysis_1_dds$Condition, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}
ddssva <- analysis_1_dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + Condition
analysis_1_dds <- DESeq(ddssva)
resultsNames(analysis_1_dds)
colSums(assay(analysis_1_dds))
colData(analysis_1_dds)
normalised_analysis_1 <- counts(analysis_1_dds, normalized=TRUE) 
head(normalised_analysis_1)
analysis_1_DEG <-results(analysis_1_dds)
analysis_1_DEGs <- analysis_1_DEG[order(analysis_1_DEG$padj),]
analysis_1_result <- as.data.frame(dplyr::mutate(as.data.frame(analysis_1_DEGs ), sig=ifelse(analysis_1_DEGs$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(analysis_1_DEGs))
sum(analysis_1_result$padj<0.05, na.rm = T)
summary(analysis_1_result)
sum(analysis_1_result$padj<0.05, na.rm = T)
print("UP:")
sum(analysis_1_result$padj<0.05 & analysis_1_result$log2FoldChange > 0, na.rm = T)
print("DOWN:")
sum(analysis_1_result$padj<0.05 & analysis_1_result$log2FoldChange < 0, na.rm = T)

# PCA
analysis_1_PCA <- vst(ddssva, blind=FALSE)
pca_1 <- plotPCA(analysis_1_PCA, intgroup=c("Condition"), returnData=TRUE) 
percentVar_1 <- round(100 * attr(pca_1, "percentVar"))

# Batch and name
ggplot(pca_1, aes(PC1, PC2, color=Condition)) + 
  geom_text_repel(aes(label=name), box.padding = 1.0, point.padding = 0.5, size = 5, max.overlaps = Inf) + 
  theme(text = element_text(size = 15)) +
  geom_point(size = 5) + ggtitle("PCA") +
  xlab(paste0("PC1: ",percentVar_1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_1[2],"% variance")) + 
  coord_fixed(ratio=1)

# Condition
ggplot(pca_1, aes(PC1, PC2, color=Condition)) + 
  theme(text = element_text(size = 15)) +
  geom_point(size = 5) + ggtitle("PCA") +
  xlab(paste0("PC1: ",percentVar_1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_1[2],"% variance")) + 
  coord_fixed(ratio=1)

# Enrichment analysis
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(ensembldb)
library(ReactomePA)
DEGs_analysis_1 <- subset(analysis_1_result, padj < 0.05)
dim(DEGs_analysis_1)
ensembl_analysis_1 <- rownames(DEGs_analysis_1)
ensembl_analysis_1 <- sub("*\\..*", "", ensembl_analysis_1)
ensembl_analysis_1
id_analysis_1 <- bitr(ensembl_analysis_1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
Enrichment_DEGs_analysis_1 <- enrichGO(gene = as.vector(id_analysis_1$ENTREZID), 
                                       keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", 
                                       pvalueCutoff=0.05, qvalueCutoff=0.05, minGSSize=2, readable = TRUE)
cluster_summary_analysis_1 <- data.frame(Enrichment_DEGs_analysis_1)
dim(cluster_summary_analysis_1)
ego_analysis_1 <- setReadable(Enrichment_DEGs_analysis_1, OrgDb = org.Hs.eg.db)
head(ego_analysis_1)
clusterProfiler::dotplot(ego_analysis_1, showCategory=20)
cnetplot(ego_analysis_1)

clusterProfiler::dotplot(ego_analysis_1, showCategory=10, font.size=23) +
  scale_y_discrete(labels = function(y) lapply(strwrap(y, width = 30, simplify = FALSE),  paste, collapse="\n")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

write.csv(ego_analysis_1, file="male_enrichment.csv")

# KEGG
kegg_organism <- "hsa"
kegg_analysis_1 <- bitr(ensembl_analysis_1, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
kegg_analysis_1 <- kegg_analysis_1[!duplicated(kegg_analysis_1[c("ENSEMBL")]),]
DEGs_analysis_1$ensembl <- rownames(DEGs_analysis_1)
DEGs_analysis_1$ensembl <- sub("*\\..*", "", DEGs_analysis_1$ensembl)
row.names(DEGs_analysis_1) <- DEGs_analysis_1$ensembl
DEGs_analysis_1$ensembl <- NULL
kegg_analysis_1_a <- DEGs_analysis_1[rownames(DEGs_analysis_1) %in% kegg_analysis_1$ENSEMBL,]
kegg_analysis_1_a$Y <- kegg_analysis_1$ENTREZID
FC_analysis_1 <- kegg_analysis_1_a$log2FoldChange
names(FC_analysis_1) <- kegg_analysis_1_a$Y
kegg_gene_list_analysis_1 <- sort(FC_analysis_1, decreasing = TRUE)
kegg_analysis_1 <- gseKEGG(geneList = kegg_gene_list_analysis_1,
                           organism     = kegg_organism,
                           nPerm        = 1000,
                           minGSSize    = 2,
                           maxGSSize    = 800,
                           pvalueCutoff = 0.10,
                           pAdjustMethod = "none",
                           keyType       = "ncbi-geneid")
gseaKEGG_results_analysis_1 <- kegg_analysis_1@result
head(gseaKEGG_results_analysis_1)
kegg_readable_analysis_1 <- setReadable(kegg_analysis_1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_cluster_summary_analysis_1 <- data.frame(kegg_readable_analysis_1)
clusterProfiler::dotplot(kegg_readable_analysis_1, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
write.csv(kegg_cluster_summary_analysis_1, file="male_KEGG.csv")

# REACTOME
reactome_analysis_1 <- enrichPathway(gene=kegg_analysis_1_a$Y, pvalueCutoff=0.05, pAdjustMethod="none", readable=T, organism="human")
reactome_analaysis_1_A <- (as.data.frame(reactome_analysis_1))
head(reactome_analysis_1)
dim(reactome_analysis_1)
clusterProfiler::dotplot(reactome_analysis_1)

clusterProfiler::dotplot(reactome_analysis_1, showCategory=10, font.size=23) +
  scale_y_discrete(labels = function(y) lapply(strwrap(y, width = 30, simplify = FALSE),  paste, collapse="\n")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

write.csv(reactome_analysis_1, file="male_reactome.csv")

# Organise data
analysis_1_merge <- merge(analysis_1_result, normalised_analysis_1, by=0, all=TRUE)
rownames(analysis_1_merge) = analysis_1_merge$Row.names
analysis_1_merge$Row.names <- NULL
head(analysis_1_merge)
analysis_1_merge_1 <- merge(analysis_1_merge, ids, by=0, all=TRUE)
rownames(analysis_1_merge_1) <- analysis_1_merge_1$Row.names
analysis_1_merge_1$Row.names <- NULL
analysis_1_merge_2 <- analysis_1_merge_1 %>% arrange(padj)
analysis_1_merge_3 <- na.omit(analysis_1_merge_2)
final_1 <- as.data.frame(dplyr::mutate(as.data.frame(analysis_1_merge_3), 
                                       M_HF =ifelse(analysis_1_merge_3$log2FoldChange>0, "UP", "DOWN")), 
                         row.names=rownames(analysis_1_merge_3))
write.csv(final_1, file="male_DEGS.csv")

# Volcano plot 
volcano_plot_1 <- analysis_1_merge_3[ which(analysis_1_merge_3$padj < 0.05),]
rownames(volcano_plot_1) = volcano_plot_1$gene_symbol

# Duplicate names  
rownames(volcano_plot_1) <- make.names(volcano_plot_1$gene_symbol, unique = TRUE) 

p <- ggplot(analysis_1_merge_3, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col = sig)) + 
  scale_color_manual(values = c("red", "black")) + 
  geom_vline(xintercept=c(-0,0), linetype="dotted") + 
  theme(legend.position = "none", text = element_text(size=20))
p + geom_text_repel(data=volcano_plot_1[1:100, ], aes(label=rownames(volcano_plot_1[1:100, ])))

##############################
# Analysis 2 (Female_HF VERSUS Female_CON) #
##############################
sample_analysis_1 <- overall_sample_conditions[which(overall_sample_conditions$Condition == "female_HF" | overall_sample_conditions$Condition == "female_CON"),]
head(sample_analysis_1)
dim(sample_analysis_1)
table(sample_analysis_1$Condition)
sample_analysis_1_ids <- sample_analysis_1$Sample_File
sample_analysis_1_ids <- as.vector(t(sample_analysis_1_ids))
head(sample_analysis_1_ids)
sample_analysis_1_count <- overall_counts[, sample_analysis_1_ids]
head(sample_analysis_1_count)
dim(sample_analysis_1_count)

### Correct batch effect using ComBat-seq
#sample_analysis_1_count <- ComBat_seq(counts = sample_analysis_1_count, batch = sample_analysis_1$Batch)

# Identify DEGs
sample_analysis_1$Batch <- factor(sample_analysis_1$Batch)
sample_analysis_1_DE <- DESeqDataSetFromMatrix(countData=sample_analysis_1_count, colData=sample_analysis_1, design=~Condition)
as.data.frame(colData(sample_analysis_1_DE))
analysis_1 <- sample_analysis_1_DE[rowSums(DESeq2::counts(sample_analysis_1_DE)) > 10, ]
analysis_1_dds <- DESeq(analysis_1)
resultsNames(analysis_1_dds)
colSums(assay(analysis_1_dds))
colData(analysis_1_dds)
normalised_analysis_1 <- counts(analysis_1_dds, normalized=TRUE) 
head(normalised_analysis_1)
analysis_1_DEG <-results(analysis_1_dds)
analysis_1_DEGs <- analysis_1_DEG[order(analysis_1_DEG$padj),]
analysis_1_result <- as.data.frame(dplyr::mutate(as.data.frame(analysis_1_DEGs ), sig=ifelse(analysis_1_DEGs$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(analysis_1_DEGs))
sum(analysis_1_result$padj<0.05, na.rm = T)
summary(analysis_1_result)
sum(analysis_1_result$padj<0.05, na.rm = T)
print("UP:")
sum(analysis_1_result$padj<0.05 & analysis_1_result$log2FoldChange > 0, na.rm = T)
print("DOWN:")
sum(analysis_1_result$padj<0.05 & analysis_1_result$log2FoldChange < 0, na.rm = T)

# Correct batch effect using sva
dat  <- counts(analysis_1_dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ Condition, colData(analysis_1_dds))
mod0 <- model.matrix(~   1, colData(analysis_1_dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
svseq$sv
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ analysis_1_dds$Condition, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}
ddssva <- analysis_1_dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + Condition
analysis_1_dds <- DESeq(ddssva)
resultsNames(analysis_1_dds)
colSums(assay(analysis_1_dds))
colData(analysis_1_dds)
normalised_analysis_1 <- counts(analysis_1_dds, normalized=TRUE) 
head(normalised_analysis_1)
analysis_1_DEG <-results(analysis_1_dds)
analysis_1_DEGs <- analysis_1_DEG[order(analysis_1_DEG$padj),]
analysis_1_result <- as.data.frame(dplyr::mutate(as.data.frame(analysis_1_DEGs ), sig=ifelse(analysis_1_DEGs$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(analysis_1_DEGs))
sum(analysis_1_result$padj<0.05, na.rm = T)
summary(analysis_1_result)
sum(analysis_1_result$padj<0.05, na.rm = T)
print("UP:")
sum(analysis_1_result$padj<0.05 & analysis_1_result$log2FoldChange > 0, na.rm = T)
print("DOWN:")
sum(analysis_1_result$padj<0.05 & analysis_1_result$log2FoldChange < 0, na.rm = T)

# PCA
analysis_1_PCA <- vst(ddssva, blind=FALSE)
pca_1 <- plotPCA(analysis_1_PCA, intgroup=c("Condition"), returnData=TRUE) 
percentVar_1 <- round(100 * attr(pca_1, "percentVar"))

# Batch and name
ggplot(pca_1, aes(PC1, PC2, color=Condition)) + 
  geom_text_repel(aes(label=name), box.padding = 1.0, point.padding = 0.5, size = 5, max.overlaps = Inf) + 
  theme(text = element_text(size = 15)) +
  geom_point(size = 5) + ggtitle("PCA") +
  xlab(paste0("PC1: ",percentVar_1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_1[2],"% variance")) + 
  coord_fixed(ratio=1)

# Condition
ggplot(pca_1, aes(PC1, PC2, color=Condition)) + 
  theme(text = element_text(size = 15)) +
  geom_point(size = 5) + ggtitle("PCA") +
  xlab(paste0("PC1: ",percentVar_1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_1[2],"% variance")) + 
  coord_fixed(ratio=1)

# Enrichment analysis
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(ensembldb)
library(ReactomePA)
DEGs_analysis_1 <- subset(analysis_1_result, padj < 0.05)
dim(DEGs_analysis_1)
ensembl_analysis_1 <- rownames(DEGs_analysis_1)
ensembl_analysis_1 <- sub("*\\..*", "", ensembl_analysis_1)
ensembl_analysis_1
id_analysis_1 <- bitr(ensembl_analysis_1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
Enrichment_DEGs_analysis_1 <- enrichGO(gene = as.vector(id_analysis_1$ENTREZID), 
                                       keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", 
                                       pvalueCutoff=0.05, qvalueCutoff=0.05, minGSSize=2, readable = TRUE)
cluster_summary_analysis_1 <- data.frame(Enrichment_DEGs_analysis_1)
dim(cluster_summary_analysis_1)
ego_analysis_1 <- setReadable(Enrichment_DEGs_analysis_1, OrgDb = org.Hs.eg.db)
head(ego_analysis_1)
clusterProfiler::dotplot(ego_analysis_1, showCategory=20)
cnetplot(ego_analysis_1)

clusterProfiler::dotplot(ego_analysis_1, showCategory=10, font.size=23) +
  scale_y_discrete(labels = function(y) lapply(strwrap(y, width = 30, simplify = FALSE),  paste, collapse="\n")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

write.csv(ego_analysis_1, file="female_enrichment.csv")

# KEGG
kegg_organism <- "hsa"
kegg_analysis_1 <- bitr(ensembl_analysis_1, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
kegg_analysis_1 <- kegg_analysis_1[!duplicated(kegg_analysis_1[c("ENSEMBL")]),]
DEGs_analysis_1$ensembl <- rownames(DEGs_analysis_1)
DEGs_analysis_1$ensembl <- sub("*\\..*", "", DEGs_analysis_1$ensembl)
row.names(DEGs_analysis_1) <- DEGs_analysis_1$ensembl
DEGs_analysis_1$ensembl <- NULL
kegg_analysis_1_a <- DEGs_analysis_1[rownames(DEGs_analysis_1) %in% kegg_analysis_1$ENSEMBL,]
kegg_analysis_1_a$Y <- kegg_analysis_1$ENTREZID
FC_analysis_1 <- kegg_analysis_1_a$log2FoldChange
names(FC_analysis_1) <- kegg_analysis_1_a$Y
kegg_gene_list_analysis_1 <- sort(FC_analysis_1, decreasing = TRUE)
kegg_analysis_1 <- gseKEGG(geneList = kegg_gene_list_analysis_1,
                           organism     = kegg_organism,
                           nPerm        = 1000,
                           minGSSize    = 2,
                           maxGSSize    = 800,
                           pvalueCutoff = 0.10,
                           pAdjustMethod = "none",
                           keyType       = "ncbi-geneid")
gseaKEGG_results_analysis_1 <- kegg_analysis_1@result
head(gseaKEGG_results_analysis_1)
kegg_readable_analysis_1 <- setReadable(kegg_analysis_1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_cluster_summary_analysis_1 <- data.frame(kegg_readable_analysis_1)
clusterProfiler::dotplot(kegg_readable_analysis_1, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
write.csv(kegg_cluster_summary_analysis_1, file="female_KEGG.csv")

# REACTOME
reactome_analysis_1 <- enrichPathway(gene=kegg_analysis_1_a$Y, pvalueCutoff=0.05, pAdjustMethod="none", readable=T, organism="human")
reactome_analaysis_1_A <- (as.data.frame(reactome_analysis_1))
head(reactome_analysis_1)
dim(reactome_analysis_1)
clusterProfiler::dotplot(reactome_analysis_1)

clusterProfiler::dotplot(reactome_analysis_1, showCategory=10, font.size=23) +
  scale_y_discrete(labels = function(y) lapply(strwrap(y, width = 30, simplify = FALSE),  paste, collapse="\n")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

write.csv(reactome_analysis_1, file="female_reactome.csv")

# Organise data
analysis_1_merge <- merge(analysis_1_result, normalised_analysis_1, by=0, all=TRUE)
rownames(analysis_1_merge) = analysis_1_merge$Row.names
analysis_1_merge$Row.names <- NULL
head(analysis_1_merge)
analysis_1_merge_1 <- merge(analysis_1_merge, ids, by=0, all=TRUE)
rownames(analysis_1_merge_1) <- analysis_1_merge_1$Row.names
analysis_1_merge_1$Row.names <- NULL
analysis_1_merge_2 <- analysis_1_merge_1 %>% arrange(padj)
analysis_1_merge_3 <- na.omit(analysis_1_merge_2)
final_1 <- as.data.frame(dplyr::mutate(as.data.frame(analysis_1_merge_3), 
                                       F_HF =ifelse(analysis_1_merge_3$log2FoldChange>0, "UP", "DOWN")), 
                         row.names=rownames(analysis_1_merge_3))
write.csv(final_1, file="female_DEGS.csv")

# Volcano plot 
volcano_plot_1 <- analysis_1_merge_3[ which(analysis_1_merge_3$padj < 0.05),]
rownames(volcano_plot_1) = volcano_plot_1$gene_symbol

# Duplicate names  
rownames(volcano_plot_1) <- make.names(volcano_plot_1$gene_symbol, unique = TRUE) 

p <- ggplot(analysis_1_merge_3, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col = sig)) + 
  scale_color_manual(values = c("red", "black")) + 
  geom_vline(xintercept=c(-0,0), linetype="dotted") + 
  theme(legend.position = "none", text = element_text(size=20))
p + geom_text_repel(data=volcano_plot_1[1:100, ], aes(label=rownames(volcano_plot_1[1:100, ])))
