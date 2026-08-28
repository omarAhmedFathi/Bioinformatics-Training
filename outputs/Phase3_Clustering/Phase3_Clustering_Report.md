# Phase 3: Identify Sub-population Structure (Clustering on PCs)

## Objective
Identify latent subpopulation structure directly from the Principal Components derived in Phase 2.

## Methodology
We applied **k-means clustering** to the top Principal Components. 

### Choosing the Number of Clusters
To determine the appropriate number of clusters ($k$), we computed the total within-cluster sum of squares (WSS) for $k=1$ to $10$ (using the elbow method). We observed the largest drop in WSS between $k=3$ and $k=4$, with the curve flattening out significantly afterwards. Therefore, we chose **$k=4$** as the optimal number of clusters for this dataset. This suggests the presence of 4 main subpopulation groups.

## Results & Comparison: 2 PCs vs. 3 PCs
We ran k-means clustering twice to compare the effect of using different numbers of PCs:
1. Using the top 2 PCs (PC1, PC2).
2. Using the top 3 PCs (PC1, PC2, PC3).

### Deliverable Plots
- **`Cluster_Plot_2PCs.pdf`**: Shows the scatter plot of PC1 vs. PC2, colored by the 4 cluster assignments.
- **`Cluster_Plot_3PCs.pdf`**: Shows a pairwise grid of scatter plots (PC1 vs. PC2, PC1 vs. PC3, and PC2 vs. PC3) using the assignments from the 3-PC clustering model.

### Comparison Note
When comparing the usage of 2 PCs versus 3 PCs for clustering, the subpopulation assignments remain highly consistent. A cross-tabulation of the cluster labels shows that the vast majority of samples cluster together identically in both approaches. Only 4 out of 156 samples (approx. 2.5%) shifted cluster memberships when the 3rd Principal Component was introduced. This indicates that the first two Principal Components capture the primary latent population structure, but the 3rd PC provides a very slight refinement in delineating the boundaries of one of the major clusters.
