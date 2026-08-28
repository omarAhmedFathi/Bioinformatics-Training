# Task 1: Population Genetics & GWAS Pipeline

This tutorial covers the first 7 phases of the bioinformatics pipeline. We start with raw genotype data, perform quality control, identify population structure, and conduct Genome-Wide Association Studies (GWAS) for both quantitative traits and binary traits. Finally, we annotate the significant findings.

---

## Phase 1: Quality Control (QC)

> [!NOTE]
> **The Biology:** Before any genomic analysis, we must ensure the data is reliable. We filter out SNPs with very low Minor Allele Frequencies (MAF), as they lack statistical power and may represent genotyping errors. We also filter out SNPs and individuals with high missingness.

### Method & Code
We use PLINK to compute allele frequencies and missingness, then visualize the distributions in R.

```R
library(ggplot2)

# Read MAF data
freq <- read.table("intermediate/phase1_freq.frq", header=TRUE, stringsAsFactors=FALSE)

cat("Min MAF:", min(freq$MAF, na.rm=TRUE), "\n")
cat("Max MAF:", max(freq$MAF, na.rm=TRUE), "\n")

# Plot MAF
ggplot(freq, aes(x=MAF)) + 
  geom_histogram(binwidth=0.01, fill="steelblue", color="black") +
  theme_minimal() +
  labs(title="Minor Allele Frequency (MAF) Distribution", x="MAF", y="Frequency")
```

### Result
Below is the distribution of Minor Allele Frequencies across our dataset. Notice how filtering will remove the left-most peak of very rare variants.

![MAF Histogram](outputs/Phase1_QC/Phase1_MAF_histogram-1.png)

---

## Phase 2 & 3: Population Structure (PCA & Clustering)

> [!TIP]
> **The Biology:** Human populations have underlying genetic structures (ancestry). If we don't account for this, our GWAS will produce false positives due to "population stratification." PCA helps us capture these ancestral axes.

### Method & Code
We calculate the eigenvectors (PCs) using PLINK, then plot them. We apply clustering (Hierarchical and DBSCAN) on the top PCs to empirically define subpopulations.

```R
pca_data <- read.table("intermediate/phase2_pca.eigenvec", header=FALSE, stringsAsFactors=FALSE)
colnames(pca_data) <- c("FID", "IID", paste0("PC", 1:10))

# Scatter plot PC1 vs PC2
ggplot(pca_data, aes(x=PC1, y=PC2)) +
  geom_point(alpha=0.7, color="darkorange") +
  theme_minimal() +
  labs(title="PCA: PC1 vs PC2", x="Principal Component 1", y="Principal Component 2")

# Hierarchical Clustering using 3 PCs
d3 <- dist(pca_data[, c("PC1", "PC2", "PC3")])
hc3 <- hclust(d3, method = "ward.D2")
pca_data$Cluster_3PC <- as.factor(cutree(hc3, k = 4))
```

### Result
The scatter plot of PC1 vs PC2 reveals distinct clusters corresponding to ancestral backgrounds. 

![PCA Plot](outputs/Phase2_PCA/Phase2_PCA_Scatter_Plot-1.png)
![Clustering Plot](outputs/Phase3_Clustering/Phase3_Cluster_Plot_3PCs-1.png)

---

## Phase 4 & 5: GWAS (Quantitative & Binary Traits)

> [!IMPORTANT]
> **The Biology:** We aim to find statistical associations between individual SNPs and a phenotype. In Phase 4, we use PC1 as a quantitative pseudo-phenotype to see which SNPs drive population structure. In Phase 5, we use Sex as a binary trait as a biological sanity check (expecting strong hits on the X/Y chromosomes).

### Method & Code
PLINK runs the linear/logistic regressions. We prepare the phenotype and covariate files in R to pass into PLINK.

```R
# Prepare phenotype for PC1
pheno_pc1 <- pca_data[, c("FID", "IID", "PC1")]
write.table(pheno_pc1, "intermediate/pheno_PC1.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Prepare covariate file (PC3 to PC10)
covar_data <- pca_data[, c("FID", "IID", paste0("PC", 3:10))]
write.table(covar_data, "intermediate/covar.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
```

*(Association testing is then executed via bash scripts using `plink --linear` and `plink --logistic`)*

---

## Phase 6 & 7: Annotation & Enrichment

> [!NOTE]
> **The Biology:** A statistically significant SNP is just a coordinate. We must annotate it (map it to a gene) and perform pathway enrichment to understand what biological pathways are actually being impacted.

### Method & Code
We use the `biomaRt` package to query Ensembl and find the nearest genes, then use `clusterProfiler` for GO (Gene Ontology) enrichment. Finally, we use `qqman` to plot beautiful Manhattan plots.

```R
library(biomaRt)
library(qqman)

# Load GWAS Results
gwas_pc1 <- read.table("outputs/Phase4_GWAS_Quantitative/Phase4_GWAS_PC1_results.txt", header=TRUE)
gwas_pc1 <- na.omit(gwas_pc1)

# Annotate SNPs using biomaRt (code abbreviated)
snp_mart <- useEnsembl(biomart="snps", dataset="hsapiens_snp")
# ... retrieval logic ...

# Draw Manhattan Plot
manhattan(gwas_pc1_annotated, chr="CHR", bp="BP", snp="Label", p="P", 
          annotatePval = 5e-8, annotateTop = FALSE, 
          suggestiveline = FALSE, main="Manhattan Plot: PC1 Association")
```

### Result
The Manhattan plot highlights the genomic regions with genome-wide significance (peaks above the threshold line).

![Manhattan PC1](outputs/Phase6_Annotation/Phase6_Manhattan_PC1_Annotated-1.png)
*(If no SNPs reached significance, a standard plot is generated).*

**Pathway Enrichment:**
Using `clusterProfiler`, we can identify which biological pathways (like immune response, metabolic regulation) these genes collectively govern.

![GO Enrichment](outputs/Phase7_Enrichment/Phase7_Enrichment_Dotplot-1.png)

---
*End of Task 1.* [Proceed to Task 2](Task2_Metabolite_GWAS.md)
