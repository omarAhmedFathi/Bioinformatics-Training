# ==========================================
# Task 2 - Phase 14: Partial Correlation Network of Metabolites
# ==========================================

cat("=== Phase 14: Partial Correlation Network Analysis ===\n\n")

if (!require("ppcor", quietly=TRUE)) install.packages("ppcor")
if (!require("igraph", quietly=TRUE)) install.packages("igraph")
if (!require("dplyr", quietly=TRUE)) install.packages("dplyr")
if (!require("ggplot2", quietly=TRUE)) install.packages("ggplot2")

library(ppcor)
library(igraph)
library(dplyr)
library(ggplot2)

# 1. Load Data
cat("Loading QC'd metabolomics dataset...\n")
qc_data <- read.csv("../intermediate/Phase11_Metabolites_QC_Final.csv", stringsAsFactors = FALSE)

# Extract only metabolite columns (drop IDs and phenotype)
metabolite_cols <- names(qc_data)[!(names(qc_data) %in% c("main_id", "mapped_id", "Diabetes", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Cluster"))]
metabolites <- qc_data[, metabolite_cols]

cat(sprintf("Number of metabolites for network: %d\n", ncol(metabolites)))

# 2. Compute Partial Correlation Matrix
cat("Computing partial correlations using ppcor...\n")
pcor_results <- pcor(metabolites)

# Extract partial correlations and p-values
pcor_est <- pcor_results$estimate
pcor_pval <- pcor_results$p.value

rownames(pcor_est) <- colnames(pcor_est) <- metabolite_cols
rownames(pcor_pval) <- colnames(pcor_pval) <- metabolite_cols

# 3. Thresholding to define edges
cat("Thresholding edges (p < 0.05 and |partial r| > 0.1)...\n")
p_val_thresh <- 0.05
r_thresh <- 0.1

# Create an adjacency matrix (1 for edge, 0 for no edge)
adj_matrix <- (pcor_pval < p_val_thresh) & (abs(pcor_est) > r_thresh)
diag(adj_matrix) <- FALSE # No self-loops

# Create the edge list
edges <- data.frame(Source = character(), Target = character(), Weight = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

for (i in 1:(ncol(pcor_est)-1)) {
  for (j in (i+1):ncol(pcor_est)) {
    if (adj_matrix[i, j]) {
      edges <- rbind(edges, data.frame(
        Source = metabolite_cols[i],
        Target = metabolite_cols[j],
        Weight = pcor_est[i, j],
        P_Value = pcor_pval[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat(sprintf("Generated %d significant edges.\n", nrow(edges)))

# 4. Build and Visualize Network Graph with igraph
cat("Building igraph network...\n")
g <- graph_from_data_frame(d = edges, directed = FALSE)

# Assign edge attributes (color based on sign of correlation, width based on strength)
E(g)$color <- ifelse(E(g)$Weight > 0, "firebrick", "steelblue")
E(g)$width <- abs(E(g)$Weight) * 10

# Calculate some network metrics for node attributes
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g)

# Save the plot
pdf("../outputs/Phase14_Network/Phase14_Partial_Correlation_Network.pdf", width = 12, height = 12)
set.seed(42)
plot(g, 
     layout = layout_with_fr(g), 
     vertex.size = 5 + (V(g)$degree * 1.5), 
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     vertex.color = "gold",
     vertex.frame.color = "gray30",
     edge.color = E(g)$color,
     edge.width = E(g)$width,
     main = "Metabolite Partial Correlation Network")
legend("topright", legend=c("Positive Correlation", "Negative Correlation"), col=c("firebrick", "steelblue"), lwd=2)
invisible(dev.off())

# 5. Export for Cytoscape
cat("Exporting node and edge tables for Cytoscape...\n")
# Node table
nodes_df <- data.frame(
  Id = V(g)$name,
  Degree = V(g)$degree,
  Betweenness = V(g)$betweenness,
  stringsAsFactors = FALSE
)

write.csv(nodes_df, "../outputs/Phase14_Network/Cytoscape_Nodes.csv", row.names = FALSE)
write.csv(edges, "../outputs/Phase14_Network/Cytoscape_Edges.csv", row.names = FALSE)
write.csv(pcor_est, "../outputs/Phase14_Network/Partial_Correlation_Matrix.csv", row.names = TRUE)

# 6. Summary Note
summary_note <- paste0(
  "Phase 14: Partial Correlation Network Analysis Summary\n\n",
  "Methodology:\n",
  "- Computed partial correlations using `ppcor` to identify direct metabolic associations while controlling for all other metabolites.\n",
  "- Thresholded edges at unadjusted p < 0.05 and |partial r| > 0.1.\n",
  "- Resulting network has ", vcount(g), " nodes (metabolites) and ", ecount(g), " edges.\n\n",
  "Top Hubs (Highest Degree):\n"
)
top_hubs <- nodes_df %>% arrange(desc(Degree)) %>% head(5)
for(i in 1:nrow(top_hubs)) {
  summary_note <- paste0(summary_note, "- ", top_hubs$Id[i], " (Degree: ", top_hubs$Degree[i], ")\n")
}
writeLines(summary_note, "../outputs/Phase14_Network/Phase14_Network_Interpretation.txt")

cat("\nPhase 14 complete. Outputs saved in ../outputs/Phase14_Network/:\n")
cat("  - Phase14_Partial_Correlation_Network.pdf\n")
cat("  - Cytoscape_Nodes.csv\n")
cat("  - Cytoscape_Edges.csv\n")
cat("  - Partial_Correlation_Matrix.csv\n")
cat("  - Phase14_Network_Interpretation.txt\n")
