# ==========================================
# Phase 1: Quality Control Plots
# ==========================================
library(ggplot2)

# Read MAF data
freq <- read.table("../intermediate/phase1_freq.frq", header=TRUE, stringsAsFactors=FALSE)

cat("Min MAF:", min(freq$MAF, na.rm=TRUE), "\n")
cat("Max MAF:", max(freq$MAF, na.rm=TRUE), "\n")

pdf("../outputs/Phase1_QC/Phase1_MAF_histogram.pdf")
print(ggplot(freq, aes(x=MAF)) + 
  geom_histogram(binwidth=0.01, fill="steelblue", color="black") +
  theme_minimal() +
  labs(title="Minor Allele Frequency (MAF) Distribution", x="MAF", y="Frequency"))
invisible(dev.off())

# Read missingness data (SNP)
lmiss <- read.table("../intermediate/phase1_missing.lmiss", header=TRUE, stringsAsFactors=FALSE)
pdf("../outputs/Phase1_QC/Phase1_SNP_missingness_histogram.pdf")
print(ggplot(lmiss, aes(x=F_MISS)) + 
  geom_histogram(binwidth=0.001, fill="darkgreen", color="black") +
  theme_minimal() +
  labs(title="Missingness per SNP", x="Fraction Missing", y="Count"))
invisible(dev.off())

# Read missingness data (Sample)
imiss <- read.table("../intermediate/phase1_missing.imiss", header=TRUE, stringsAsFactors=FALSE)
pdf("../outputs/Phase1_QC/Phase1_Sample_missingness_histogram.pdf")
print(ggplot(imiss, aes(x=F_MISS)) + 
  geom_histogram(binwidth=0.001, fill="darkred", color="black") +
  theme_minimal() +
  labs(title="Missingness per Sample", x="Fraction Missing", y="Count"))
invisible(dev.off())

# ==========================================
# Phase 2: PCA Plots
# ==========================================

# Read Eigenvectors for Scatter Plot
pca_data <- read.table("../intermediate/phase2_pca.eigenvec", header=FALSE, stringsAsFactors=FALSE)
colnames(pca_data) <- c("FID", "IID", paste0("PC", 1:10))

# Scatter plot PC1 vs PC2
pdf("../outputs/Phase2_PCA/Phase2_PCA_Scatter_Plot.pdf")
print(ggplot(pca_data, aes(x=PC1, y=PC2)) +
  geom_point(alpha=0.7, color="darkorange") +
  theme_minimal() +
  labs(title="PCA: PC1 vs PC2", x="Principal Component 1", y="Principal Component 2"))
invisible(dev.off())

# Read Eigenvalues for Scree Plot
eigenvalues <- read.table("../intermediate/phase2_pca.eigenval", header=FALSE)$V1
pve <- eigenvalues / sum(eigenvalues) * 100 # Note: This is variance explained relative to top 10 PCs

scree_data <- data.frame(
  PC = factor(1:length(eigenvalues)),
  Eigenvalue = eigenvalues,
  PVE = pve
)

# Scree plot
pdf("../outputs/Phase2_PCA/Phase2_PCA_Scree_Plot.pdf")
print(ggplot(scree_data, aes(x=PC, y=Eigenvalue)) +
  geom_bar(stat="identity", fill="steelblue", color="black") +
  geom_line(aes(x=as.numeric(PC), y=Eigenvalue), color="red", linewidth=1) +
  geom_point(aes(x=as.numeric(PC), y=Eigenvalue), color="red", size=2) +
  theme_minimal() +
  labs(title="Scree Plot (Top 10 PCs)", x="Principal Component", y="Eigenvalue (Variance)"))
invisible(dev.off())

# ==========================================
# Phase 3: Clustering on Principal Components
# ==========================================
library(gridExtra)

# Hierarchical Clustering using 2 PCs
d2 <- dist(pca_data[, c("PC1", "PC2")])
hc2 <- hclust(d2, method = "ward.D2")
pca_data$Cluster_2PC <- as.factor(cutree(hc2, k = 4))

pdf("../outputs/Phase3_Clustering/Phase3_Cluster_Plot_2PCs.pdf", width=6, height=5)
print(ggplot(pca_data, aes(x=PC1, y=PC2, color=Cluster_2PC)) +
  geom_point(alpha=0.7, size=2) +
  theme_minimal() +
  labs(title="Hierarchical Clustering (k=4) on 2 PCs", color="Cluster"))
invisible(dev.off())

# Hierarchical Clustering using 3 PCs
d3 <- dist(pca_data[, c("PC1", "PC2", "PC3")])
hc3 <- hclust(d3, method = "ward.D2")
pca_data$Cluster_3PC <- as.factor(cutree(hc3, k = 4))

# Plot Dendrogram
pdf("../outputs/Phase3_Clustering/Phase3_Dendrogram.pdf", width=10, height=5)
plot(hc3, labels=FALSE, main="Dendrogram (Ward's method on 3 PCs)", xlab="", sub="")
rect.hclust(hc3, k=4, border="red")
invisible(dev.off())

p1 <- ggplot(pca_data, aes(x=PC1, y=PC2, color=Cluster_3PC)) + geom_point(alpha=0.7) + theme_minimal() + theme(legend.position="none")
p2 <- ggplot(pca_data, aes(x=PC1, y=PC3, color=Cluster_3PC)) + geom_point(alpha=0.7) + theme_minimal() + theme(legend.position="none")
p3 <- ggplot(pca_data, aes(x=PC2, y=PC3, color=Cluster_3PC)) + geom_point(alpha=0.7) + theme_minimal()

pdf("../outputs/Phase3_Clustering/Phase3_Cluster_Plot_3PCs.pdf", width=12, height=4)
grid.arrange(p1, p2, p3, ncol=3, top="Hierarchical Clustering (k=4) on 3 PCs")
invisible(dev.off())

# DBSCAN Clustering using 3 PCs
library(dbscan)
kNN_dists <- kNNdist(pca_data[, c("PC1", "PC2", "PC3")], k = 5)
eps_val <- quantile(kNN_dists, 0.9) # Dynamically choose epsilon
set.seed(42)
db <- dbscan(pca_data[, c("PC1", "PC2", "PC3")], eps = eps_val, minPts = 5)
pca_data$Cluster_DBSCAN <- as.factor(db$cluster)

p1_db <- ggplot(pca_data, aes(x=PC1, y=PC2, color=Cluster_DBSCAN)) + geom_point(alpha=0.7) + theme_minimal() + theme(legend.position="none")
p2_db <- ggplot(pca_data, aes(x=PC1, y=PC3, color=Cluster_DBSCAN)) + geom_point(alpha=0.7) + theme_minimal() + theme(legend.position="none")
p3_db <- ggplot(pca_data, aes(x=PC2, y=PC3, color=Cluster_DBSCAN)) + geom_point(alpha=0.7) + theme_minimal()

pdf("../outputs/Phase3_Clustering/Phase3_DBSCAN_Plot_3PCs.pdf", width=12, height=4)
grid.arrange(p1_db, p2_db, p3_db, ncol=3, top="DBSCAN Clustering on 3 PCs")
invisible(dev.off())

# Save cross-tabulation
write.table(table(pca_data$Cluster_2PC, pca_data$Cluster_3PC), "../intermediate/cluster_compare.txt", quote=FALSE)

# ==========================================
# Phase 4: Prepare GWAS (Quantitative Trait)
# ==========================================

# Phenotype files
pheno_pc1 <- pca_data[, c("FID", "IID", "PC1")]
pheno_pc2 <- pca_data[, c("FID", "IID", "PC2")]

write.table(pheno_pc1, "../intermediate/pheno_PC1.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(pheno_pc2, "../intermediate/pheno_PC2.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Covariate file (PC3 to PC10)
covar_data <- pca_data[, c("FID", "IID", paste0("PC", 3:10))]
write.table(covar_data, "../intermediate/covar.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)

# ==========================================
# Phase 5: Prepare GWAS (Sex)
# ==========================================

fam_data <- read.table("../data/Qatari156_filtered_pruned.fam", header=FALSE, stringsAsFactors=FALSE)

cat("Sex coding counts:\n")
print(table(fam_data$V5))

pheno_sex <- data.frame(FID=fam_data$V1, IID=fam_data$V2, Pheno=fam_data$V5)
write.table(pheno_sex, "../intermediate/pheno_sex.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

# ==========================================
# Phase 6: Gene and functional Annotation of Significant SNPs
# ==========================================
library(biomaRt)
library(qqman)
library(dplyr)

# 1. Load GWAS Results
gwas_pc1 <- read.table("../outputs/Phase4_GWAS_Quantitative/Phase4_GWAS_PC1_results.txt", header=TRUE, stringsAsFactors=FALSE)
gwas_sex <- read.table("../outputs/Phase5_GWAS_Sex/Phase5_GWAS_Sex_results.txt", header=TRUE, stringsAsFactors=FALSE)

# Ensure P-values are numeric and omit NAs
gwas_pc1$P <- as.numeric(gwas_pc1$P)
gwas_sex$P <- as.numeric(gwas_sex$P)
gwas_pc1 <- na.omit(gwas_pc1)
gwas_sex <- na.omit(gwas_sex)

# 2. Identify Significant SNPs (P < 5e-8)
sig_pc1 <- gwas_pc1[gwas_pc1$P < 5e-8, ]
sig_sex <- gwas_sex[gwas_sex$P < 5e-8, ]
all_sig_snps <- unique(c(sig_pc1$SNP, sig_sex$SNP))

if (length(all_sig_snps) > 0) {
  # 3. Annotate SNPs using biomaRt
  snp_mart <- useEnsembl(biomart="snps", dataset="hsapiens_snp")
  annotations <- getBM(
    attributes = c("refsnp_id", "ensembl_gene_stable_id", "consequence_type_tv", "distance_to_transcript"), 
    filters = "snp_filter", 
    values = all_sig_snps, 
    mart = snp_mart
  )
  
  # Get Gene Symbols from Ensembl
  gene_mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = unique(annotations$ensembl_gene_stable_id),
    mart = gene_mart
  )
  
  # Merge annotations with gene symbols
  annotated_snps <- merge(annotations, gene_info, by.x="ensembl_gene_stable_id", by.y="ensembl_gene_id", all.x=TRUE)
  write.csv(annotated_snps, "../outputs/Phase6_Annotation/Phase6_Annotated_Significant_SNPs.csv", row.names=FALSE)
  
  # 4. Gene-labeled Manhattan Plots
  label_data <- annotated_snps %>%
    group_by(refsnp_id) %>%
    summarise(Gene = first(na.omit(hgnc_symbol))) %>%
    filter(Gene != "")
  
  # Merge with GWAS results to replace SNP ID with Gene for significant SNPs
  gwas_pc1_annotated <- merge(gwas_pc1, label_data, by.x="SNP", by.y="refsnp_id", all.x=TRUE)
  gwas_pc1_annotated$Label <- ifelse(is.na(gwas_pc1_annotated$Gene), as.character(gwas_pc1_annotated$SNP), as.character(gwas_pc1_annotated$Gene))
  
  gwas_sex_annotated <- merge(gwas_sex, label_data, by.x="SNP", by.y="refsnp_id", all.x=TRUE)
  gwas_sex_annotated$Label <- ifelse(is.na(gwas_sex_annotated$Gene), as.character(gwas_sex_annotated$SNP), as.character(gwas_sex_annotated$Gene))
  
  # Manhattan Plot for PC1
  pdf("../outputs/Phase6_Annotation/Phase6_Manhattan_PC1_Annotated.pdf", width=10, height=6)
  manhattan(gwas_pc1_annotated, chr="CHR", bp="BP", snp="Label", p="P", 
            annotatePval = 5e-8, annotateTop = FALSE, 
            suggestiveline = FALSE, main="Manhattan Plot: PC1 Association")
  invisible(dev.off())
  
  # Manhattan Plot for Sex
  pdf("../outputs/Phase6_Annotation/Phase6_Manhattan_Sex_Annotated.pdf", width=10, height=6)
  manhattan(gwas_sex_annotated, chr="CHR", bp="BP", snp="Label", p="P", 
            annotatePval = 5e-8, annotateTop = FALSE, 
            suggestiveline = FALSE, main="Manhattan Plot: Sex Association")
  invisible(dev.off())
} else {
  cat("No significant SNPs found with P < 5e-8. Generating standard Manhattan plots without annotations.\n")
  pdf("../outputs/Phase6_Annotation/Phase6_Manhattan_PC1.pdf", width=10, height=6)
  manhattan(gwas_pc1, suggestiveline = FALSE, main="Manhattan Plot: PC1 Association")
  invisible(dev.off())
  
  pdf("../outputs/Phase6_Annotation/Phase6_Manhattan_Sex.pdf", width=10, height=6)
  manhattan(gwas_sex, suggestiveline = FALSE, main="Manhattan Plot: Sex Association")
  invisible(dev.off())
}

# ==========================================
# Phase 7: Pathway Enrichment Analysis
# ==========================================
library(clusterProfiler)
library(org.Hs.eg.db)

if (exists("all_sig_snps") && length(all_sig_snps) > 0 && exists("label_data") && nrow(label_data) > 0) {
  genes <- unique(label_data$Gene)
  
  # Convert HGNC symbols to Entrez IDs (required for clusterProfiler KEGG/GO)
  gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_entrez) > 0) {
    # Run GO Enrichment (Biological Process)
    ego <- enrichGO(gene          = gene_entrez$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    
    # Save results table
    if (!is.null(ego) && nrow(ego) > 0) {
      write.csv(as.data.frame(ego), "../outputs/Phase7_Enrichment/Phase7_Enriched_GO_Pathways.csv", row.names=FALSE)
      
      # Plotting
      pdf("../outputs/Phase7_Enrichment/Phase7_Enrichment_Dotplot.pdf", width=10, height=8)
      print(dotplot(ego, showCategory=15, title="GO Enrichment Analysis (BP)"))
      invisible(dev.off())
      
      pdf("../outputs/Phase7_Enrichment/Phase7_Enrichment_Barplot.pdf", width=10, height=8)
      print(barplot(ego, showCategory=15, title="GO Enrichment Analysis (BP)"))
      invisible(dev.off())
    } else {
      cat("No significantly enriched GO pathways found.\n")
    }
  } else {
    cat("Could not map Gene Symbols to Entrez IDs.\n")
  }
} else {
  cat("No significant genes available for pathway enrichment analysis.\n")
}
