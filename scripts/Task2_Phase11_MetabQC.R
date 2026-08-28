# ==========================================
# Task 2 - Phase 11: Metabolite Data Quality Control
# ==========================================

cat("=== Phase 11: Metabolite Data Quality Control ===\n\n")

# Load required packages
if (!require("dplyr", quietly=TRUE)) install.packages("dplyr")
library(dplyr)

# 1. Load Data
metab_raw <- read.csv("../data/Qatari_metabolomics(in).csv", stringsAsFactors = FALSE)
mapping <- read.csv("../data/mapping(Sheet1).csv", stringsAsFactors = FALSE)

# Merge ID mapping to get the main_id (QBC-xxx) 
# so it matches Task 1's IDs
metab <- merge(mapping, metab_raw, by = "mapped_id")

# Save a copy of the diabetes status and IDs
sample_info <- metab[, c("main_id", "mapped_id", "Diabetes")]

# Extract only the metabolite columns (drop IDs and Diabetes)
metabolites <- metab[, !(colnames(metab) %in% c("main_id", "mapped_id", "Diabetes"))]

cat(sprintf("Initial dataset: %d samples, %d metabolites.\n\n", nrow(metabolites), ncol(metabolites)))

# 2. Metabolite-level missingness (>20%)
metab_missingness <- colSums(is.na(metabolites)) / nrow(metabolites)
metab_miss_table <- data.frame(Metabolite = names(metab_missingness), Missing_Prop = metab_missingness)

metabs_to_remove_missing <- names(metab_missingness[metab_missingness > 0.20])
metabolites_qc1 <- metabolites[, !(names(metabolites) %in% metabs_to_remove_missing)]

cat(sprintf("Step 1: Removed %d metabolites with >20%% missingness.\n", length(metabs_to_remove_missing)))

# 3. Sample-level missingness (>20%)
sample_missingness <- rowSums(is.na(metabolites_qc1)) / ncol(metabolites_qc1)
sample_miss_table <- data.frame(Sample = sample_info$main_id, Missing_Prop = sample_missingness)

samples_to_keep <- sample_missingness <= 0.20
metabolites_qc2 <- metabolites_qc1[samples_to_keep, ]
sample_info_qc <- sample_info[samples_to_keep, ]

cat(sprintf("Step 2: Removed %d samples with >20%% missingness.\n", sum(!samples_to_keep)))

# 4. Remove uninformative metabolites (zero or near-zero variance)
# Calculate variance for each metabolite (removing NAs)
variances <- sapply(metabolites_qc2, var, na.rm = TRUE)
# Let's consider var < 1e-6 as near zero variance
metabs_to_remove_var <- names(variances[variances < 1e-6 | is.na(variances)])
metabolites_qc3 <- metabolites_qc2[, !(names(metabolites_qc2) %in% metabs_to_remove_var)]

cat(sprintf("Step 3: Removed %d metabolites with near-zero variance.\n", length(metabs_to_remove_var)))

# Total removed metabolites list
all_removed_metabs <- unique(c(metabs_to_remove_missing, metabs_to_remove_var))
passed_metabs <- names(metabolites_qc3)

# 5. Normalize/standardize metabolite values (z-score scaling)
metabolites_final <- as.data.frame(scale(metabolites_qc3))

# Re-attach sample info
final_dataset <- cbind(sample_info_qc, metabolites_final)

cat(sprintf("\nFinal QC'd dataset: %d samples, %d metabolites.\n", nrow(metabolites_final), ncol(metabolites_final)))

# 6. Save deliverables
write.csv(metab_miss_table, "../intermediate/Phase11_Metabolite_Missingness.csv", row.names = FALSE)
write.csv(sample_miss_table, "../intermediate/Phase11_Sample_Missingness.csv", row.names = FALSE)
write.table(all_removed_metabs, "../intermediate/Phase11_Removed_Metabolites.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(passed_metabs, "../intermediate/Phase11_QC_Passed_Metabolites.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.csv(final_dataset, "../intermediate/Phase11_Metabolites_QC_Final.csv", row.names = FALSE)

cat("\nPhase 11 complete. Outputs saved in intermediate_files/:\n")
cat("  - Phase11_Metabolite_Missingness.csv\n")
cat("  - Phase11_Sample_Missingness.csv\n")
cat("  - Phase11_Removed_Metabolites.txt\n")
cat("  - Phase11_QC_Passed_Metabolites.txt\n")
cat("  - Phase11_Metabolites_QC_Final.csv (ready for Phase 12)\n")
