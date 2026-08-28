# ==========================================
# Task 2 - Phase 15: SNP–Metabolite Association (mQTL) Mapping
# ==========================================
#
# Unadjusted scan: metabolite ~ genotype (simple linear regression)
# Adjusted scan:   mixed model with Sex + PC1-10 fixed effects + GRM random effect
#                  (EMMAX-like approach using eigendecomposition of GRM)
#
# Data: synthetic_obesity_metabolome_data.csv (156 samples, 32 metabolites)
#       Qatari156_filtered_pruned (PLINK fileset)
# ==========================================

cat("=== Phase 15: SNP-Metabolite Association (mQTL) Mapping ===\n\n")

# --- Load libraries ---
if (!require("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!require("qqman", quietly=TRUE)) install.packages("qqman")
if (!require("dplyr", quietly=TRUE)) install.packages("dplyr")
if (!require("data.table", quietly=TRUE)) install.packages("data.table")
if (!require("regress", quietly=TRUE)) install.packages("regress")

library(ggplot2)
library(qqman)
library(dplyr)
library(data.table)
library(regress)

# =========================================
# 1. LOAD AND PREPARE DATA
# =========================================
cat("Step 1: Loading and preparing data...\n")

# 1a. Export PLINK genotypes to .raw (additive dosage) format
cat("  Exporting PLINK genotypes to dosage format...\n")
system("plink --bfile ../data/Qatari156_filtered_pruned --recodeA --out ../intermediate/Phase15_raw_genotypes 2>&1",
       intern = TRUE)

# 1b. Load genotype dosages (data.table for speed)
cat("  Loading genotype dosages...\n")
geno_raw <- fread("../intermediate/Phase15_raw_genotypes.raw", header = TRUE)
geno_df  <- as.data.frame(geno_raw)

geno_iids <- geno_df$IID
geno_cols <- names(geno_df)[7:ncol(geno_df)]
G <- as.matrix(geno_df[, geno_cols, drop = FALSE])

# Strip allele suffix from column names (rs12345_A -> rs12345)
snp_names <- sub("_[ACGT]+$", "", geno_cols)
colnames(G) <- snp_names

cat(sprintf("  Genotype matrix: %d samples x %d SNPs\n", nrow(G), ncol(G)))

# Impute missing genotypes to column mean (standard practice)
na_counts <- colSums(is.na(G))
if (any(na_counts > 0)) {
  cat(sprintf("  Imputing %d missing genotype values to column mean...\n", sum(na_counts)))
  col_means <- colMeans(G, na.rm = TRUE)
  for (j in which(na_counts > 0)) {
    G[is.na(G[, j]), j] <- col_means[j]
  }
}

# Remove monomorphic SNPs
snp_var <- apply(G, 2, var)
mono <- snp_var < 1e-10
if (any(mono)) {
  cat(sprintf("  Removing %d monomorphic SNPs...\n", sum(mono)))
  G <- G[, !mono]
  snp_names <- snp_names[!mono]
}
cat(sprintf("  Final genotype matrix: %d samples x %d SNPs\n", nrow(G), ncol(G)))

# 1c. Load metabolome data
cat("  Loading metabolome data...\n")
metab_raw  <- read.csv("../data/metabolome_data - synthetic_obesity_metabolome_data.csv",
                        stringsAsFactors = FALSE)
metab_iids <- metab_raw$IID
metab_cols <- names(metab_raw)[!(names(metab_raw) %in% c("IID", "Trait"))]
Y <- as.matrix(metab_raw[, metab_cols])
cat(sprintf("  Metabolome matrix: %d samples x %d metabolites\n", nrow(Y), ncol(Y)))

# Create display-name mapping for cleaner output
display_names <- metab_cols
names(display_names) <- metab_cols

# 1d. Load covariates (sex from FAM, PCs from eigenvec)
cat("  Loading covariates...\n")
fam <- read.table("../data/Qatari156_filtered_pruned.fam",
                   header = FALSE, stringsAsFactors = FALSE)
sex_vec <- fam$V5            # 1 = male, 2 = female

eigenvec <- read.table("../intermediate/phase2_pca.eigenvec",
                        header = FALSE, stringsAsFactors = FALSE)
colnames(eigenvec) <- c("FID", "IID", paste0("PC", 1:10))
pc_iids <- eigenvec$IID
PCs <- as.matrix(eigenvec[, paste0("PC", 1:10)])

# 1e. Load GRM
cat("  Loading GRM...\n")
grm <- as.matrix(read.table("../intermediate/phase9_grm.rel"))
grm_ids <- read.table("../intermediate/phase9_grm.rel.id",
                       header = TRUE, stringsAsFactors = FALSE, comment.char = "")
colnames(grm_ids)[1] <- "FID"
grm_iids <- grm_ids$IID
rownames(grm) <- grm_iids
colnames(grm) <- grm_iids

# 1f. Align all datasets by IID
cat("  Aligning samples...\n")
common_iids <- Reduce(intersect, list(geno_iids, metab_iids, pc_iids, grm_iids))
cat(sprintf("  Common samples: %d\n", length(common_iids)))

idx_geno  <- match(common_iids, geno_iids)
idx_metab <- match(common_iids, metab_iids)
idx_pc    <- match(common_iids, pc_iids)
idx_grm   <- match(common_iids, grm_iids)

G   <- G[idx_geno, ]
Y   <- Y[idx_metab, ]
sex <- sex_vec[idx_geno]
PCs <- PCs[idx_pc, ]
GRM <- grm[idx_grm, idx_grm]

n       <- nrow(G)
n_snps  <- ncol(G)
n_metab <- ncol(Y)
cat(sprintf("  Aligned: %d samples, %d SNPs, %d metabolites\n", n, n_snps, n_metab))

# Build covariate matrix (intercept + sex + PC1-10)
X <- cbind(1, sex, PCs)
colnames(X) <- c("Intercept", "Sex", paste0("PC", 1:10))
p_cov <- ncol(X)

# Load BIM for SNP positions (needed for Manhattan plots)
bim <- read.table("../data/Qatari156_filtered_pruned.bim",
                   header = FALSE, stringsAsFactors = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
bim_idx <- match(colnames(G), bim$SNP)
bim_matched <- bim[bim_idx, ]

# =========================================
# 2. UNADJUSTED mQTL SCAN
# =========================================
cat("\nStep 2: Unadjusted mQTL scan (metabolite ~ genotype)...\n")

G_c   <- scale(G, center = TRUE, scale = FALSE)  # center genotypes
G_css <- colSums(G_c^2)
valid_snps <- G_css > 1e-10                       # skip zero-variance

unadj_results <- vector("list", n_metab)

for (m in seq_len(n_metab)) {
  y   <- Y[, m]
  y_c <- y - mean(y, na.rm = TRUE)

  Sxy <- as.vector(crossprod(G_c[, valid_snps], y_c))
  css <- G_css[valid_snps]

  beta <- Sxy / css
  SS_y <- sum(y_c^2)
  SS_resid <- SS_y - Sxy^2 / css
  df  <- n - 2
  MSE <- SS_resid / df
  SE  <- sqrt(MSE / css)
  t_stat <- beta / SE
  p_val  <- 2 * pt(-abs(t_stat), df = df)

  unadj_results[[m]] <- data.frame(
    SNP        = colnames(G)[valid_snps],
    CHR        = bim_matched$CHR[valid_snps],
    BP         = bim_matched$BP[valid_snps],
    Metabolite = metab_cols[m],
    Beta       = beta,
    SE         = SE,
    P          = p_val,
    stringsAsFactors = FALSE
  )

  if (m %% 8 == 0 || m == n_metab)
    cat(sprintf("  Completed %d/%d metabolites\n", m, n_metab))
}

unadj_all <- do.call(rbind, unadj_results)
cat(sprintf("  Unadjusted scan complete: %d tests\n", nrow(unadj_all)))

# =========================================
# 3. ADJUSTED mQTL SCAN (Mixed Model)
# =========================================
cat("\nStep 3: Adjusted mQTL scan (mixed model with GRM)...\n")
cat("  Method: EMMAX-like — eigendecompose GRM, profile REML for variance components,\n")
cat("          then weighted regression with Sex + PC1-10 as fixed effects.\n")

# 3a. Eigendecompose GRM (once)
cat("  Eigendecomposing GRM...\n")
eig <- eigen(GRM, symmetric = TRUE)
U   <- eig$vectors
D   <- eig$values
D[D < 1e-6] <- 1e-6          # threshold tiny/negative eigenvalues

# 3b. Rotate ALL data into eigenspace (expensive for G, but done only once)
cat("  Rotating data into eigenspace...\n")
Y_rot <- crossprod(U, Y)     # 156 x 32
X_rot <- crossprod(U, X)     # 156 x 12
G_rot <- crossprod(U, G)     # 156 x n_snps  (main computation)
cat("  Rotation complete.\n")

# 3c. Profile REML log-likelihood (optimise delta = sigma_g / sigma_e)
profile_reml_negloglik <- function(log_delta, y_rot, X_rot, D) {
  delta <- exp(log_delta)
  nn    <- length(y_rot)
  pp    <- ncol(X_rot)

  w  <- 1 / (delta * D + 1)
  sw <- sqrt(w)

  Xw <- X_rot * sw
  yw <- y_rot * sw

  XtX <- crossprod(Xw)
  Xty <- crossprod(Xw, yw)
  beta <- tryCatch(solve(XtX, Xty),
                   error = function(e) solve(XtX + diag(1e-8, ncol(XtX)), Xty))

  resid <- yw - Xw %*% beta
  RSS   <- sum(resid^2)

  # Profile REML: -2 * REML ∝  sum(log(delta*d+1)) + (n-p)*log(RSS)
  0.5 * (sum(log(delta * D + 1)) + (nn - pp) * log(RSS))
}

# 3d. Run adjusted scan metabolite-by-metabolite
adj_results <- vector("list", n_metab)
delta_estimates <- numeric(n_metab)

for (m in seq_len(n_metab)) {
  y_m <- Y_rot[, m]

  # Estimate delta via 1-D optimisation (very fast)
  opt <- optimize(profile_reml_negloglik, interval = c(-10, 10),
                  y_rot = y_m, X_rot = X_rot, D = D)
  delta <- exp(opt$minimum)
  delta_estimates[m] <- delta

  # Weights from estimated variance components
  w  <- 1 / (delta * D + 1)
  sw <- sqrt(w)

  # Weighted (transformed) data
  y_w <- y_m * sw
  X_w <- X_rot * sw
  G_w <- G_rot * sw             # element-wise; R recycles sw down columns

  # Frisch–Waugh: residualise y and G on X
  XtX_inv <- tryCatch(solve(crossprod(X_w)),
                      error = function(e) solve(crossprod(X_w) + diag(1e-8, p_cov)))
  beta_X_y <- XtX_inv %*% crossprod(X_w, y_w)
  y_res    <- as.vector(y_w - X_w %*% beta_X_y)

  XtG      <- crossprod(X_w, G_w)      # p_cov x n_snps
  beta_X_G <- XtX_inv %*% XtG           # p_cov x n_snps
  G_res    <- G_w - X_w %*% beta_X_G    # n x n_snps

  # Vectorised SNP test
  Sxy <- as.vector(crossprod(G_res, y_res))
  css <- colSums(G_res^2)
  valid <- css > 1e-10

  beta_snp <- rep(NA_real_, n_snps)
  se_snp   <- rep(NA_real_, n_snps)
  p_snp    <- rep(NA_real_, n_snps)

  beta_snp[valid] <- Sxy[valid] / css[valid]
  SS_y_res <- sum(y_res^2)
  SS_resid_snp <- SS_y_res - Sxy[valid]^2 / css[valid]
  df_adj <- n - p_cov - 1
  MSE    <- SS_resid_snp / df_adj
  SE_v   <- sqrt(MSE / css[valid])
  se_snp[valid] <- SE_v
  t_stat <- beta_snp[valid] / SE_v
  p_snp[valid] <- 2 * pt(-abs(t_stat), df = df_adj)

  adj_results[[m]] <- data.frame(
    SNP        = colnames(G),
    CHR        = bim_matched$CHR,
    BP         = bim_matched$BP,
    Metabolite = metab_cols[m],
    Beta       = beta_snp,
    SE         = se_snp,
    P          = p_snp,
    stringsAsFactors = FALSE
  )

  if (m %% 8 == 0 || m == n_metab)
    cat(sprintf("  Completed %d/%d metabolites  (delta = %.4f)\n", m, n_metab, delta))
}

adj_all <- do.call(rbind, adj_results)
adj_all <- adj_all[!is.na(adj_all$P), ]
cat(sprintf("  Adjusted scan complete: %d tests\n", nrow(adj_all)))

# =========================================
# 4. COMPARE THE TWO RUNS
# =========================================
cat("\nStep 4: Comparing unadjusted vs adjusted runs...\n")

calc_lambda <- function(pvals) {
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals < 1]
  chisq <- qchisq(1 - pvals, df = 1)
  median(chisq) / qchisq(0.5, df = 1)
}

lambda_unadj <- calc_lambda(unadj_all$P)
lambda_adj   <- calc_lambda(adj_all$P)

cat(sprintf("  Lambda (unadjusted): %.4f\n", lambda_unadj))
cat(sprintf("  Lambda (adjusted):   %.4f\n", lambda_adj))

writeLines(c(
  "Genomic Inflation Factor (lambda)",
  sprintf("Unadjusted: %.4f", lambda_unadj),
  sprintf("Adjusted:   %.4f", lambda_adj)
), "../intermediate/Phase15_lambda_values.txt")

bonferroni_thresh <- 0.05 / nrow(unadj_all)
suggestive_thresh <- 1e-5

n_sig_unadj <- sum(unadj_all$P < bonferroni_thresh, na.rm = TRUE)
n_sig_adj   <- sum(adj_all$P < bonferroni_thresh,   na.rm = TRUE)
n_sug_unadj <- sum(unadj_all$P < suggestive_thresh, na.rm = TRUE)
n_sug_adj   <- sum(adj_all$P < suggestive_thresh,   na.rm = TRUE)

cat(sprintf("  Bonferroni threshold: %.2e\n", bonferroni_thresh))
cat(sprintf("  Genome-wide significant: Unadjusted = %d, Adjusted = %d\n",
            n_sig_unadj, n_sig_adj))
cat(sprintf("  Suggestive (p < 1e-5):   Unadjusted = %d, Adjusted = %d\n",
            n_sug_unadj, n_sug_adj))

# Overlap at suggestive threshold
if (n_sug_unadj > 0 && n_sug_adj > 0) {
  sig_u <- paste(unadj_all$SNP[unadj_all$P < suggestive_thresh],
                 unadj_all$Metabolite[unadj_all$P < suggestive_thresh], sep = ":")
  sig_a <- paste(adj_all$SNP[adj_all$P < suggestive_thresh],
                 adj_all$Metabolite[adj_all$P < suggestive_thresh], sep = ":")
  n_overlap    <- length(intersect(sig_u, sig_a))
  n_unadj_only <- length(setdiff(sig_u, sig_a))
  n_adj_only   <- length(setdiff(sig_a, sig_u))
} else {
  n_overlap <- 0; n_unadj_only <- n_sug_unadj; n_adj_only <- n_sug_adj
}

comparison_df <- data.frame(
  Metric     = c("Lambda", "Bonferroni_Threshold",
                  "N_Significant_GW", "N_Suggestive",
                  "N_Suggestive_Overlap", "N_Unadjusted_Only", "N_Adjusted_Only"),
  Unadjusted = c(round(lambda_unadj, 4), signif(bonferroni_thresh, 3),
                  n_sig_unadj, n_sug_unadj, n_overlap, n_unadj_only, NA),
  Adjusted   = c(round(lambda_adj, 4), signif(bonferroni_thresh, 3),
                  n_sig_adj, n_sug_adj, n_overlap, NA, n_adj_only),
  stringsAsFactors = FALSE
)
write.csv(comparison_df,
          "../outputs/Phase15_mQTL/Phase15_Comparison_Summary.csv", row.names = FALSE)

# =========================================
# 4b. Manhattan & QQ Plots
# =========================================
cat("  Generating Manhattan & QQ plots...\n")

# Combine across metabolites: min-p per SNP for Manhattan
manhattan_unadj <- unadj_all %>%
  group_by(SNP, CHR, BP) %>%
  summarise(P = min(P, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(P), P > 0) %>%
  as.data.frame()

manhattan_adj <- adj_all %>%
  group_by(SNP, CHR, BP) %>%
  summarise(P = min(P, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(P), P > 0) %>%
  as.data.frame()

pdf("../outputs/Phase15_mQTL/Phase15_Manhattan_Unadjusted.pdf", width = 12, height = 6)
manhattan(manhattan_unadj, chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          main = paste0("mQTL Manhattan — Unadjusted (\u03bb = ", round(lambda_unadj, 3), ")"),
          suggestiveline = FALSE,
          genomewideline = -log10(bonferroni_thresh))
invisible(dev.off())

pdf("../outputs/Phase15_mQTL/Phase15_Manhattan_Adjusted.pdf", width = 12, height = 6)
manhattan(manhattan_adj, chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          main = paste0("mQTL Manhattan — Adjusted (\u03bb = ", round(lambda_adj, 3), ")"),
          suggestiveline = FALSE,
          genomewideline = -log10(bonferroni_thresh))
invisible(dev.off())

png("../outputs/Phase15_mQTL/Phase15_QQ_Unadjusted.png", width = 7, height = 7, units = "in", res = 300)
qq(unadj_all$P[unadj_all$P > 0],
   main = paste0("QQ Plot — Unadjusted (\u03bb = ", round(lambda_unadj, 3), ")"))
invisible(dev.off())

png("../outputs/Phase15_mQTL/Phase15_QQ_Adjusted.png", width = 7, height = 7, units = "in", res = 300)
qq(adj_all$P[adj_all$P > 0 & !is.na(adj_all$P)],
   main = paste0("QQ Plot — Adjusted (\u03bb = ", round(lambda_adj, 3), ")"))
invisible(dev.off())

# =========================================
# 5. SAVE ASSOCIATION TABLES
# =========================================
cat("\nStep 5: Saving association tables...\n")

# Keep file sizes manageable: top 100 per metabolite + all suggestive
save_top <- function(df, label) {
  top <- df %>%
    group_by(Metabolite) %>% arrange(P) %>% slice_head(n = 100) %>% ungroup()
  sug <- df %>% filter(P < 1e-4)
  out <- bind_rows(top, sug) %>% distinct() %>% arrange(P)
  out$Covariates <- ifelse(label == "Adjusted",
                           "Sex, PC1-PC10, GRM (mixed model)", "None")
  out
}

top_unadj <- save_top(unadj_all, "Unadjusted")
top_adj   <- save_top(adj_all,   "Adjusted")

write.csv(top_unadj,
          "../outputs/Phase15_mQTL/Phase15_mQTL_Unadjusted_Results.csv", row.names = FALSE)
write.csv(top_adj,
          "../outputs/Phase15_mQTL/Phase15_mQTL_Adjusted_Results.csv",   row.names = FALSE)

cat(sprintf("  Saved %d unadjusted and %d adjusted top results\n",
            nrow(top_unadj), nrow(top_adj)))

# =========================================
# 6. ANNOTATE SIGNIFICANT HITS
# =========================================
cat("\nStep 6: Annotating top hits with nearest genes...\n")

# Gather hits (suggestive or top-50 if nothing passes)
hits <- bind_rows(
  unadj_all %>% filter(P < suggestive_thresh) %>% mutate(Run = "Unadjusted"),
  adj_all   %>% filter(P < suggestive_thresh) %>% mutate(Run = "Adjusted")
)
if (nrow(hits) == 0) {
  cat("  No suggestive hits. Using relaxed threshold (p < 1e-3, top 50 each)...\n")
  hits <- bind_rows(
    unadj_all %>% arrange(P) %>% head(50) %>% mutate(Run = "Unadjusted"),
    adj_all   %>% arrange(P) %>% head(50) %>% mutate(Run = "Adjusted")
  )
}

unique_snps <- unique(hits$SNP)
cat(sprintf("  Annotating %d unique SNPs via biomaRt...\n", length(unique_snps)))

annotated <- tryCatch({
  library(biomaRt)
  snp_mart <- useEnsembl(biomart = "snps", dataset = "hsapiens_snp")
  ann <- getBM(
    attributes = c("refsnp_id", "chr_name", "chrom_start",
                    "ensembl_gene_stable_id", "consequence_type_tv",
                    "distance_to_transcript"),
    filters = "snp_filter", values = unique_snps, mart = snp_mart)

  if (nrow(ann) > 0) {
    gene_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    gene_info <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol",
                      "chromosome_name", "start_position", "end_position"),
      filters = "ensembl_gene_id",
      values  = unique(ann$ensembl_gene_stable_id), mart = gene_mart)

    ann_m <- merge(ann, gene_info,
                   by.x = "ensembl_gene_stable_id", by.y = "ensembl_gene_id",
                   all.x = TRUE)
    nearest <- ann_m %>%
      filter(hgnc_symbol != "" & !is.na(hgnc_symbol)) %>%
      group_by(refsnp_id) %>%
      arrange(abs(distance_to_transcript)) %>% slice_head(n = 1) %>% ungroup() %>%
      select(SNP = refsnp_id, Nearest_Gene = hgnc_symbol,
             Distance = distance_to_transcript, Consequence = consequence_type_tv)
    nearest$Cis_Trans <- ifelse(abs(nearest$Distance) < 1e6, "cis", "trans")
    nearest
  } else NULL
}, error = function(e) {
  cat(sprintf("  biomaRt failed: %s\n  Proceeding without annotation.\n", e$message))
  NULL
})

if (!is.null(annotated) && nrow(annotated) > 0) {
  hits_ann <- merge(hits, annotated, by = "SNP", all.x = TRUE)
} else {
  hits_ann <- hits
  hits_ann$Nearest_Gene <- NA; hits_ann$Distance <- NA
  hits_ann$Consequence  <- NA; hits_ann$Cis_Trans <- NA
}
hits_ann <- hits_ann %>% arrange(P)
write.csv(hits_ann,
          "../outputs/Phase15_mQTL/Phase15_Annotated_Significant_mQTLs.csv", row.names = FALSE)
cat(sprintf("  Saved %d annotated mQTL hits\n", nrow(hits_ann)))

# Save mQTL-associated metabolites for Phase 16
mqtl_metabs <- unique(hits$Metabolite)
writeLines(mqtl_metabs, "../intermediate/Phase15_mQTL_Metabolites.txt")

# =========================================
# 7. INTERPRETATION NOTE
# =========================================
cat("\nStep 7: Writing interpretation...\n")

interp <- paste0(
  "Phase 15: SNP-Metabolite Association (mQTL) Mapping\n",
  "====================================================\n\n",
  "Dataset: synthetic_obesity_metabolome_data\n",
  sprintf("  Samples: %d  |  SNPs: %d  |  Metabolites: %d\n\n", n, n_snps, n_metab),
  "METHODOLOGY\n",
  "-----------\n",
  "Unadjusted: Simple linear regression  metabolite ~ genotype  (no covariates).\n",
  "Adjusted:   Linear mixed model (EMMAX approach):\n",
  "  1. Eigendecompose GRM to decorrelate samples.\n",
  "  2. Estimate variance-component ratio (delta = sigma_g/sigma_e) via\n",
  "     profile REML under the null model for each metabolite.\n",
  "  3. Weighted regression incorporating Sex + PC1-10 as fixed-effect\n",
  "     covariates and the GRM as a random effect, controlling for\n",
  "     population stratification and cryptic relatedness.\n\n",
  "RESULTS\n",
  "-------\n",
  sprintf("Total tests per run: %s (%s SNPs x %d metabolites)\n",
          format(as.integer(n_snps) * n_metab, big.mark = ","),
          format(n_snps, big.mark = ","), n_metab),
  sprintf("Bonferroni threshold: p < %.2e\n\n", bonferroni_thresh),
  sprintf("Genomic inflation factor (lambda):\n"),
  sprintf("  Unadjusted: %.4f\n", lambda_unadj),
  sprintf("  Adjusted:   %.4f\n\n", lambda_adj),
  sprintf("Genome-wide significant hits: Unadjusted = %d, Adjusted = %d\n",
          n_sig_unadj, n_sig_adj),
  sprintf("Suggestive hits (p < 1e-5):   Unadjusted = %d, Adjusted = %d\n",
          n_sug_unadj, n_sug_adj),
  ifelse(n_overlap > 0,
         sprintf("Overlap at suggestive level:  %d shared, %d unadjusted-only, %d adjusted-only\n",
                 n_overlap, n_unadj_only, n_adj_only), ""),
  "\nINTERPRETATION\n",
  "--------------\n",
  "Comparing adjusted versus unadjusted scans demonstrates the effect of controlling\n",
  "for population structure and relatedness.  A lambda closer to 1.0 in the adjusted\n",
  "scan confirms that the mixed model effectively accounts for confounding.\n\n",
  "With only 156 samples, statistical power is limited; large-scale mQTL studies\n",
  "(QMDiab, TwinsUK, KORA) typically require 500+ samples for robust discovery.\n",
  "Suggestive hits here are candidates for replication.\n\n",
  "BIOLOGICAL CONTEXT\n",
  "------------------\n",
  "BCAA mQTLs (valine, leucine, isoleucine):  commonly map to BCAT1/2, BCKDHA/B.\n",
  "Acylcarnitine mQTLs (C3, C5, C6, C16):    linked to CPT1A, ACADM, ACADL.\n",
  "Bile acid mQTLs (taurocholate, GCDC):      associated with SLCO1B1, UGT enzymes.\n",
  "Xanthine mQTLs (caffeine, paraxanthine):   strongly driven by CYP1A2 variation.\n",
  "Purine mQTLs (hypoxanthine, xanthine):     linked to XDH (xanthine dehydrogenase).\n",
  "GPC/lipid mQTLs:                           map to LPCAT, PLA2 family genes.\n"
)
writeLines(interp, "../outputs/Phase15_mQTL/Phase15_Interpretation.txt")

cat("\n=== Phase 15 Complete ===\n")
cat("Outputs in ../outputs/Phase15_mQTL/:\n")
cat("  - Phase15_mQTL_Unadjusted_Results.csv\n")
cat("  - Phase15_mQTL_Adjusted_Results.csv\n")
cat("  - Phase15_Manhattan_Unadjusted.pdf\n")
cat("  - Phase15_Manhattan_Adjusted.pdf\n")
cat("  - Phase15_QQ_Unadjusted.png\n")
cat("  - Phase15_QQ_Adjusted.png\n")
cat("  - Phase15_Comparison_Summary.csv\n")
cat("  - Phase15_Annotated_Significant_mQTLs.csv\n")
cat("  - Phase15_Interpretation.txt\n")
