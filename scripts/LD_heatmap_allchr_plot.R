# ==========================================
# Per-Chromosome LD Heatmaps (all chromosomes)
# Uses raster-based image() for large matrices
# ==========================================

library(ggplot2)

indir  <- "../intermediate/ld_decay/per_chr"
outdir <- "../outputs/Phase8_LD"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

chromosomes <- c(1:22, 23)
chr_labels  <- c(as.character(1:22), "X")

# Color palette (white -> yellow -> orange -> dark red)
ld_colors <- colorRampPalette(c("white", "#FFF7BC", "#FEE391",
                                 "#FEC44F", "#FE9929", "#D95F0E", "#993404"))(256)

# -------------------------------------------
# 1. Individual per-chromosome heatmaps (combined in one big PDF)
# -------------------------------------------
pdf(paste0(outdir, "/Phase8_LD_Heatmap_AllChromosomes.pdf"), width=14, height=12)

for (i in seq_along(chromosomes)) {
  chr <- chromosomes[i]
  lbl <- chr_labels[i]
  
  ld_file <- paste0(indir, "/chr", chr, "_ld.ld")
  snp_file <- paste0(indir, "/chr", chr, "_ld.snplist")
  
  if (!file.exists(ld_file)) {
    cat("Skipping Chr", lbl, "- file not found\n")
    next
  }
  
  cat("Plotting Chr", lbl, "...")
  
  ld_mat <- as.matrix(read.table(ld_file))
  n_snps <- nrow(ld_mat)
  
  # Cap values at 1
  ld_mat[ld_mat > 1] <- 1
  
  # Use image() for large matrices — much faster than ggplot geom_tile
  par(mar=c(4, 4, 4, 5))
  image(1:n_snps, 1:n_snps, ld_mat,
        col=ld_colors,
        zlim=c(0, 1),
        xlab="SNP index (genomic order)",
        ylab="SNP index (genomic order)",
        main=paste0("LD Heatmap — Chromosome ", lbl, " (", n_snps, " SNPs)"),
        useRaster=TRUE,
        axes=FALSE)
  
  # Add axis ticks at regular intervals
  tick_at <- pretty(1:n_snps, n=8)
  tick_at <- tick_at[tick_at >= 1 & tick_at <= n_snps]
  axis(1, at=tick_at, cex.axis=0.8)
  axis(2, at=tick_at, cex.axis=0.8)
  box()
  
  # Color bar legend
  legend_y <- seq(1, n_snps, length.out=6)
  legend_labels <- sprintf("%.1f", seq(0, 1, by=0.2))
  
  # Summary stats
  lower_vals <- ld_mat[lower.tri(ld_mat)]
  high_ld <- sum(lower_vals > 0.8)
  mod_ld  <- sum(lower_vals > 0.5 & lower_vals <= 0.8)
  
  mtext(paste0("Mean r² = ", round(mean(lower_vals), 4),
               " | Max r² = ", round(max(lower_vals), 4),
               " | Pairs r²>0.8: ", high_ld,
               " | Pairs r²>0.5: ", high_ld + mod_ld),
        side=1, line=2.8, cex=0.8, col="grey30")
  
  cat(" done (", n_snps, "SNPs, max r²=", round(max(lower_vals), 4), ")\n")
  
  rm(ld_mat, lower_vals)
  gc(verbose=FALSE)
}

invisible(dev.off())

# -------------------------------------------
# 2. Summary panel: all chromosomes side by side (smaller thumbnails)
# -------------------------------------------
cat("\n=== Creating summary panel ===\n")

pdf(paste0(outdir, "/Phase8_LD_Heatmap_Summary_Panel.pdf"), width=20, height=16)

par(mfrow=c(5, 5), mar=c(2, 2, 2.5, 1), oma=c(2, 2, 3, 0))

for (i in seq_along(chromosomes)) {
  chr <- chromosomes[i]
  lbl <- chr_labels[i]
  
  ld_file <- paste0(indir, "/chr", chr, "_ld.ld")
  
  if (!file.exists(ld_file)) next
  
  ld_mat <- as.matrix(read.table(ld_file))
  ld_mat[ld_mat > 1] <- 1
  n_snps <- nrow(ld_mat)
  
  image(1:n_snps, 1:n_snps, ld_mat,
        col=ld_colors,
        zlim=c(0, 1),
        xlab="", ylab="",
        main=paste0("Chr ", lbl, " (", n_snps, ")"),
        cex.main=1.1,
        useRaster=TRUE,
        axes=FALSE)
  box()
  
  rm(ld_mat)
  gc(verbose=FALSE)
}

# Fill remaining panels with empty plots
remaining <- 25 - length(chromosomes)
if (remaining > 0) {
  for (j in 1:remaining) plot.new()
}

mtext("Genome-wide LD Heatmaps (All Chromosomes) — Qatari Cohort, LD-Pruned Dataset",
      outer=TRUE, cex=1.3, font=2, line=1)

invisible(dev.off())

cat("\n=== All chromosome LD heatmaps complete ===\n")
cat("Outputs:\n")
cat("  1. Phase8_LD_Heatmap_AllChromosomes.pdf  (one page per chromosome, full detail)\n")
cat("  2. Phase8_LD_Heatmap_Summary_Panel.pdf   (all chromosomes on one page)\n")
