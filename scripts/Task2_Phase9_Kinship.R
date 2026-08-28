# ==========================================
# Task 2 - Phase 9: Kinship Matrix and Relative Counting
# ==========================================

library(ggplot2)
library(reshape2)
library(dplyr)

cat("=== Phase 9: Kinship Matrix & Relative Counting ===\n\n")

# -------------------------------------------
# 1. Load KING kinship table
# -------------------------------------------
king <- read.table("../intermediate/phase9_king.kin0", header=TRUE, stringsAsFactors=FALSE, comment.char="")
# Fix the #FID1 column name
colnames(king)[1] <- "FID1"
cat("Kinship table loaded:", nrow(king), "pairs\n")
cat("Columns:", paste(colnames(king), collapse=", "), "\n")
cat("Kinship range:", min(king$KINSHIP), "to", max(king$KINSHIP), "\n\n")

# -------------------------------------------
# 2. Classify pairs by KING thresholds
# -------------------------------------------
classify_kinship <- function(k) {
  case_when(
    k > 0.354   ~ "Duplicate/MZ twin",
    k > 0.177   ~ "1st-degree",
    k > 0.0884  ~ "2nd-degree",
    k > 0.0442  ~ "3rd-degree",
    TRUE         ~ "Unrelated"
  )
}

king$Relationship <- classify_kinship(king$KINSHIP)
king$Relationship <- factor(king$Relationship,
  levels=c("Duplicate/MZ twin", "1st-degree", "2nd-degree", "3rd-degree", "Unrelated"))

# -------------------------------------------
# 3. Report relatedness counts
# -------------------------------------------
rel_table <- as.data.frame(table(king$Relationship))
colnames(rel_table) <- c("Relatedness_Category", "Number_of_Pairs")

cat("=== Relatedness Classification (KING Thresholds) ===\n")
cat(sprintf("%-20s  %-15s  %s\n", "Category", "Kinship Range", "Pairs"))
cat(sprintf("%-20s  %-15s  %s\n", "--------------------", "---------------", "-----"))
cat(sprintf("%-20s  %-15s  %d\n", "Duplicate/MZ twin", "> 0.354", rel_table$Number_of_Pairs[rel_table$Relatedness_Category=="Duplicate/MZ twin"]))
cat(sprintf("%-20s  %-15s  %d\n", "1st-degree", "0.177 - 0.354", rel_table$Number_of_Pairs[rel_table$Relatedness_Category=="1st-degree"]))
cat(sprintf("%-20s  %-15s  %d\n", "2nd-degree", "0.0884 - 0.177", rel_table$Number_of_Pairs[rel_table$Relatedness_Category=="2nd-degree"]))
cat(sprintf("%-20s  %-15s  %d\n", "3rd-degree", "0.0442 - 0.0884", rel_table$Number_of_Pairs[rel_table$Relatedness_Category=="3rd-degree"]))
cat(sprintf("%-20s  %-15s  %d\n", "Unrelated", "< 0.0442", rel_table$Number_of_Pairs[rel_table$Relatedness_Category=="Unrelated"]))

# Count non-unrelated pairs and unique individuals
non_unrelated <- king[king$Relationship != "Unrelated", ]
n_relative_pairs <- nrow(non_unrelated)
individuals_in_pairs <- unique(c(non_unrelated$IID1, non_unrelated$IID2))
n_individuals <- length(individuals_in_pairs)

cat(sprintf("\n=== Summary ===\n"))
cat(sprintf("Total pairs analyzed: %d\n", nrow(king)))
cat(sprintf("Non-unrelated pairs: %d\n", n_relative_pairs))
cat(sprintf("Individuals in at least one non-unrelated pair: %d (out of 156)\n", n_individuals))

# List the most related pairs
if (n_relative_pairs > 0) {
  cat("\n=== Most Related Pairs (top 20 by kinship) ===\n")
  top_pairs <- non_unrelated %>% arrange(desc(KINSHIP)) %>% head(20)
  for (i in 1:nrow(top_pairs)) {
    cat(sprintf("  %s - %s : kinship = %.4f (%s)\n",
      top_pairs$IID1[i], top_pairs$IID2[i], 
      top_pairs$KINSHIP[i], as.character(top_pairs$Relationship[i])))
  }
}

# Save relatedness table
write.csv(rel_table, "../intermediate/Phase9_Relatedness_Table.csv", row.names=FALSE)

# -------------------------------------------
# 4. Kinship Heatmap
# -------------------------------------------

# Load the GRM for the full heatmap (more informative than KING for visualization)
grm <- as.matrix(read.table("../intermediate/phase9_grm.rel"))
grm_ids <- read.table("../intermediate/phase9_grm.rel.id", header=TRUE, stringsAsFactors=FALSE, comment.char="")
colnames(grm_ids)[1] <- "FID"
sample_ids <- grm_ids$IID
rownames(grm) <- sample_ids
colnames(grm) <- sample_ids

cat("\nGRM loaded:", nrow(grm), "x", ncol(grm), "\n")
cat("GRM diagonal range:", min(diag(grm)), "-", max(diag(grm)), "\n")
cat("GRM off-diagonal range:", min(grm[lower.tri(grm)]), "-", max(grm[lower.tri(grm)]), "\n")

# Order samples by hierarchical clustering for better heatmap structure
hc <- hclust(as.dist(1 - grm), method="ward.D2")
order_idx <- hc$order
grm_ordered <- grm[order_idx, order_idx]

# Melt for ggplot
grm_long <- melt(grm_ordered)
colnames(grm_long) <- c("Sample1", "Sample2", "Relatedness")

# Maintain clustering order
grm_long$Sample1 <- factor(grm_long$Sample1, levels=sample_ids[order_idx])
grm_long$Sample2 <- factor(grm_long$Sample2, levels=sample_ids[order_idx])

# Main deliverable: Kinship Heatmap
pdf("../outputs/Phase9_Kinship/Phase9_Kinship_Heatmap.pdf", width=14, height=12)
print(
  ggplot(grm_long, aes(x=Sample1, y=Sample2, fill=Relatedness)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B"),
      values = scales::rescale(c(min(grm_long$Relatedness), 0, 0.05, 0.1, 0.25, max(grm_long$Relatedness))),
      name = "GRM\nRelatedness"
    ) +
    theme_minimal(base_size=10) +
    theme(
      axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=3),
      axis.text.y = element_text(size=3),
      plot.title = element_text(hjust=0.5, face="bold", size=14),
      plot.subtitle = element_text(hjust=0.5, size=10),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = "Kinship Heatmap (GRM)",
      subtitle = paste0("156 Qatari Samples - Hierarchically Clustered (", n_relative_pairs, " non-unrelated pairs)"),
      x = "Sample", y = "Sample"
    ) +
    coord_fixed()
)
invisible(dev.off())

# -------------------------------------------
# 5. Distribution of kinship coefficients
# -------------------------------------------
pdf("../intermediate/Phase9_Kinship_Distribution.pdf", width=10, height=6)
print(
  ggplot(king, aes(x=KINSHIP)) +
    geom_histogram(binwidth=0.005, fill="steelblue", color="black", alpha=0.8) +
    geom_vline(xintercept=0.354, color="red", linetype="dashed", linewidth=0.8) +
    geom_vline(xintercept=0.177, color="orange", linetype="dashed", linewidth=0.8) +
    geom_vline(xintercept=0.0884, color="gold", linetype="dashed", linewidth=0.8) +
    geom_vline(xintercept=0.0442, color="green", linetype="dashed", linewidth=0.8) +
    annotate("text", x=0.354, y=Inf, label="Dup/MZ", vjust=-0.5, hjust=-0.1, color="red", size=3) +
    annotate("text", x=0.177, y=Inf, label="1st", vjust=-0.5, hjust=-0.1, color="orange", size=3) +
    annotate("text", x=0.0884, y=Inf, label="2nd", vjust=-0.5, hjust=-0.1, color="gold3", size=3) +
    annotate("text", x=0.0442, y=Inf, label="3rd", vjust=-0.5, hjust=-0.1, color="green4", size=3) +
    theme_minimal() +
    labs(
      title = "Distribution of KING Kinship Coefficients",
      subtitle = "Dashed lines show standard relatedness thresholds",
      x = "KING Kinship Coefficient", y = "Count"
    )
)
invisible(dev.off())

# -------------------------------------------
# 6. Exclusion recommendation for Phase 11
# -------------------------------------------
cat("\n=== Recommendation for Phase 11 (Mixed-Model GWAS) ===\n")
if (n_relative_pairs > 0) {
  dup_pairs <- sum(king$Relationship == "Duplicate/MZ twin")
  first_pairs <- sum(king$Relationship == "1st-degree")
  if (dup_pairs > 0) {
    cat("WARNING: Duplicate/MZ twin pairs detected. One sample from each pair\n")
    cat("should be removed before Phase 11.\n")
  }
  if (first_pairs > 0) {
    cat("NOTE: 1st-degree relatives detected. These are handled by the mixed model\n")
    cat("(GRM accounts for relatedness), but consider excluding if sample size allows.\n")
  }
  cat("\nFor Phase 11, the mixed model with GRM will account for the kinship\n")
  cat("structure. No exclusions are strictly necessary unless duplicates exist,\n")
  cat("since the GRM is designed to correct for relatedness.\n")
} else {
  cat("No closely related pairs detected. No exclusions needed for Phase 11.\n")
  cat("The GRM will still capture subtle population structure.\n")
}

cat("\nPhase 9 complete. Outputs:\n")
cat("  - Phase9_Kinship_Heatmap.pdf (main deliverable)\n")
cat("  - intermediate_files/Phase9_Kinship_Distribution.pdf (supplementary)\n")
cat("  - intermediate_files/Phase9_Relatedness_Table.csv (table)\n")
cat("  - intermediate_files/phase9_grm.rel (GRM for Phase 10 & 11)\n")
