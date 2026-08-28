# ==========================================
# Genome-wide LD Decay Plot
# Plots average r² vs physical distance across all chromosomes
# ==========================================

library(ggplot2)
library(dplyr)

cat("=== Loading LD data ===\n")
ld <- read.table("../intermediate/ld_decay/genomewide_ld.ld", header=TRUE)
cat("Total LD pairs loaded:", nrow(ld), "\n")

# Compute physical distance in bp
ld$DIST <- abs(ld$BP_B - ld$BP_A)

# Remove self-pairs (distance = 0)
ld <- ld[ld$DIST > 0, ]
cat("Pairs after removing self-hits:", nrow(ld), "\n")
cat("Distance range:", min(ld$DIST), "-", max(ld$DIST), "bp\n")
cat("Mean r²:", round(mean(ld$R2), 5), "\n")
cat("Median r²:", round(median(ld$R2), 5), "\n\n")

# -------------------------------------------
# 1. Bin by distance and compute mean r²
# -------------------------------------------
# Use variable bin sizes: fine bins for short distances, coarser for long
breaks <- c(seq(0, 50000, by=5000),
            seq(75000, 200000, by=25000),
            seq(250000, 500000, by=50000))

ld$dist_bin <- cut(ld$DIST, breaks=breaks, include.lowest=TRUE)

ld_summary <- ld %>%
  filter(!is.na(dist_bin)) %>%
  group_by(dist_bin) %>%
  summarise(
    mean_r2 = mean(R2, na.rm=TRUE),
    median_r2 = median(R2, na.rm=TRUE),
    n_pairs = n(),
    mean_dist_kb = mean(DIST) / 1000,
    .groups = "drop"
  )

cat("=== LD Decay Summary by Distance Bin ===\n")
print(as.data.frame(ld_summary), row.names=FALSE)

# -------------------------------------------
# 2. Main LD Decay Plot (mean r² vs distance)
# -------------------------------------------
dir.create("../outputs/Phase8_LD", showWarnings=FALSE, recursive=TRUE)

pdf("../outputs/Phase8_LD/Phase8_LD_Decay_Genomewide.pdf", width=12, height=7)

print(
  ggplot(ld_summary, aes(x=mean_dist_kb, y=mean_r2)) +
    geom_point(aes(size=n_pairs), color="#2196f3", alpha=0.8) +
    geom_line(color="#1565c0", linewidth=1) +
    geom_smooth(method="loess", se=TRUE, color="#d32f2f", linetype="dashed",
                fill="#ffcdd2", alpha=0.3, linewidth=0.8) +
    scale_size_continuous(name="Pair Count", labels=scales::comma) +
    scale_x_continuous(breaks=seq(0, 500, by=50)) +
    theme_minimal(base_size=13) +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=16),
      plot.subtitle = element_text(hjust=0.5, size=11, color="grey40"),
      panel.grid.minor = element_blank(),
      legend.position = c(0.85, 0.75),
      legend.background = element_rect(fill="white", color="grey80")
    ) +
    labs(
      title = "Genome-wide LD Decay Plot",
      subtitle = paste0("Qatari Cohort (n=156) | 67,735 pruned SNPs | ",
                        scales::comma(nrow(ld)), " pairwise comparisons (≤500 kb)"),
      x = "Physical Distance (kb)",
      y = expression(paste("Mean ", r^2))
    ) +
    annotate("text", x=350, y=max(ld_summary$mean_r2)*0.9,
             label=paste0("Overall mean r² = ", round(mean(ld$R2), 5)),
             hjust=0, size=4, color="grey30")
)

invisible(dev.off())

# -------------------------------------------
# 3. Per-chromosome LD decay
# -------------------------------------------
ld$CHR <- ld$CHR_A

ld_chr_summary <- ld %>%
  mutate(dist_bin_kb = floor(DIST / 10000) * 10) %>%  # 10kb bins
  group_by(CHR, dist_bin_kb) %>%
  summarise(
    mean_r2 = mean(R2, na.rm=TRUE),
    n_pairs = n(),
    .groups = "drop"
  ) %>%
  filter(n_pairs >= 10)  # require at least 10 pairs per bin

pdf("../outputs/Phase8_LD/Phase8_LD_Decay_PerChromosome.pdf", width=16, height=10)

print(
  ggplot(ld_chr_summary, aes(x=dist_bin_kb, y=mean_r2, color=factor(CHR))) +
    geom_line(alpha=0.7, linewidth=0.6) +
    facet_wrap(~paste0("Chr ", CHR), scales="free_y", ncol=5) +
    theme_minimal(base_size=10) +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=14),
      plot.subtitle = element_text(hjust=0.5, size=10, color="grey40"),
      legend.position = "none",
      strip.text = element_text(face="bold"),
      panel.grid.minor = element_blank()
    ) +
    scale_color_manual(values=rep(c("#1565c0","#d32f2f","#2e7d32","#f57f17",
                                     "#6a1b9a","#00838f","#e65100"), 4)) +
    labs(
      title = "LD Decay by Chromosome",
      subtitle = "Qatari Cohort — LD-pruned dataset (10 kb bins, ≥10 pairs per bin)",
      x = "Physical Distance (kb)",
      y = expression(paste("Mean ", r^2))
    )
)

invisible(dev.off())

# -------------------------------------------
# 4. Scatter density plot (all pairs)
# -------------------------------------------
# Sample for scatter to avoid overplotting
set.seed(42)
ld_sample <- ld[sample(nrow(ld), min(50000, nrow(ld))), ]

pdf("../outputs/Phase8_LD/Phase8_LD_Scatter_Density.pdf", width=12, height=7)

print(
  ggplot(ld_sample, aes(x=DIST/1000, y=R2)) +
    geom_point(alpha=0.05, size=0.3, color="#1565c0") +
    geom_smooth(method="loess", se=TRUE, color="#d32f2f", linewidth=1,
                fill="#ffcdd2", alpha=0.4) +
    theme_minimal(base_size=13) +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=16),
      plot.subtitle = element_text(hjust=0.5, size=11, color="grey40"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Pairwise LD vs Physical Distance (Scatter)",
      subtitle = paste0("Random sample of 50,000 SNP pairs | Qatari Cohort"),
      x = "Physical Distance (kb)",
      y = expression(r^2)
    )
)

invisible(dev.off())

cat("\n=== LD Decay Analysis Complete ===\n")
cat("Outputs:\n")
cat("  1. Phase8_LD_Decay_Genomewide.pdf  — Main LD decay curve\n")
cat("  2. Phase8_LD_Decay_PerChromosome.pdf — Per-chromosome LD decay\n")
cat("  3. Phase8_LD_Scatter_Density.pdf — Scatter density of all pairs\n")
