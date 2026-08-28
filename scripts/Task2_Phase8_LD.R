# ==========================================
# Task 2 - Phase 8: Linkage Disequilibrium (LD) Heatmap
# Generates ONLY two clean PDF heatmaps for the two regions:
#   1. Phase8_LD_Heatmap_Chr6_HLA.pdf (Reference HLA/MHC region)
#   2. Phase8_LD_Heatmap_Chr11_TopHit.pdf (Task 1 Top GWAS Hit rs10466604)
# ==========================================

library(ggplot2)
library(reshape2)

cat("=== Task 2 - Phase 8: LD Heatmap Analysis ===\n\n")

# Load BIM for positional info
bim <- read.table("../data/Qatari156_filtered_pruned.bim", header=FALSE, stringsAsFactors=FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

# Helper function to generate PLINK LD matrix if missing
generate_plink_ld <- function(chr, from_bp, to_bp, out_prefix) {
  ld_file <- paste0(out_prefix, ".ld")
  if (!file.exists(ld_file)) {
    cat(sprintf("Generating PLINK LD matrix for Chr %s (%d - %d bp)...\n", chr, from_bp, to_bp))
    cmd <- sprintf("plink --bfile ../data/Qatari156_filtered_pruned --chr %s --from-bp %d --to-bp %d --r2 square --write-snplist --out %s",
                   chr, from_bp, to_bp, out_prefix)
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  }
}

# Ensure PLINK outputs exist for both regions
generate_plink_ld(6, 25000000, 26000000, "../intermediate/phase8_ld")
generate_plink_ld(11, 123500000, 125000000, "../intermediate/phase8_top_hit_chr11_ld")

# Helper function to generate single LD Heatmap PDF
plot_ld_heatmap <- function(ld_path, snp_list_path, region_title, region_subtitle, pdf_out) {
  cat(sprintf("Generating LD Heatmap: %s -> %s\n", region_title, pdf_out))
  
  ld_matrix <- as.matrix(read.table(ld_path))
  snp_names <- read.table(snp_list_path, header=FALSE, stringsAsFactors=FALSE)$V1
  
  region_bim <- bim[bim$SNP %in% snp_names, ]
  region_bim <- region_bim[match(snp_names, region_bim$SNP), ]
  
  rownames(ld_matrix) <- snp_names
  colnames(ld_matrix) <- snp_names
  
  ld_long <- melt(ld_matrix)
  colnames(ld_long) <- c("SNP1", "SNP2", "r2")
  ld_long$SNP1 <- factor(ld_long$SNP1, levels=snp_names)
  ld_long$SNP2 <- factor(ld_long$SNP2, levels=snp_names)
  
  pdf(pdf_out, width=14, height=12)
  p <- ggplot(ld_long, aes(x=SNP1, y=SNP2, fill=r2)) +
    geom_tile(color="grey90", linewidth=0.1) +
    scale_fill_gradientn(
      colors = c("white", "#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#D95F0E", "#993404"),
      values = c(0, 0.1, 0.2, 0.3, 0.5, 0.7, 1),
      limits = c(0, 1),
      name = expression(r^2)
    ) +
    theme_minimal(base_size=10) +
    theme(
      axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=6),
      axis.text.y = element_text(size=6),
      plot.title = element_text(hjust=0.5, face="bold", size=14),
      plot.subtitle = element_text(hjust=0.5, size=10),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = paste("Linkage Disequilibrium (LD) Heatmap -", region_title),
      subtitle = region_subtitle,
      x = "SNP (ordered by genomic position)",
      y = "SNP (ordered by genomic position)"
    ) +
    coord_fixed()
  print(p)
  invisible(dev.off())
}

# Clean existing directory to leave ONLY the two requested PDFs
output_dir <- "../outputs/Phase8_LD"
old_pdfs <- list.files(output_dir, pattern = "\\.pdf$", full.names = TRUE)
if (length(old_pdfs) > 0) {
  file.remove(old_pdfs)
}

# 1. Generate PDF for Chromosome 6 HLA Region
plot_ld_heatmap(
  "../intermediate/phase8_ld.ld",
  "../intermediate/phase8_ld.snplist",
  "Chromosome 6 HLA Region",
  "Chromosome 6: 25.0-26.0 Mb - 37 SNPs",
  file.path(output_dir, "Phase8_LD_Heatmap_Chr6_HLA.pdf")
)

# 2. Generate PDF for Chromosome 11 Top Hit Region
plot_ld_heatmap(
  "../intermediate/phase8_top_hit_chr11_ld.ld",
  "../intermediate/phase8_top_hit_chr11_ld.snplist",
  "Chromosome 11 Top Hit Region",
  "Chromosome 11: 123.5-125.0 Mb (Lead SNP rs10466604, P = 2.1e-27) - 44 SNPs",
  file.path(output_dir, "Phase8_LD_Heatmap_Chr11_TopHit.pdf")
)

cat("\nPhase 8 Complete. Outputs in ../outputs/Phase8_LD/:\n")
cat("  1. Phase8_LD_Heatmap_Chr6_HLA.pdf\n")
cat("  2. Phase8_LD_Heatmap_Chr11_TopHit.pdf\n")
