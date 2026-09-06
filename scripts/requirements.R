# ==========================================
# Bioinformatics Portfolio — Package Installer
# ==========================================
# Run this script to install all required R packages.
# Usage: Rscript requirements.R
# ==========================================

cat("Installing required CRAN packages...\n")
cran_packages <- c(
  # Data manipulation
  "tidyverse", "dplyr", "data.table", "reshape2",
  # Visualization
  "ggplot2", "qqman", "pheatmap", "gridExtra", "igraph",
  # Statistics & ML
  "randomForest", "caret", "glmnet", "ppcor", "regress",
  # Clustering
  "dbscan",
  # Mixed models
  "sommer"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

cat("\nInstalling Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c(
  "biomaRt",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE)
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

cat("\n✅ All packages installed successfully!\n")
cat("\nExternal tools required (install separately):\n")
cat("  - PLINK 1.9: https://www.cog-genomics.org/plink/1.9/\n")
cat("  - GCTA: https://yanglab.westlake.edu.cn/software/gcta/\n")
cat("  - KING: https://www.kingrelatedness.com/\n")
