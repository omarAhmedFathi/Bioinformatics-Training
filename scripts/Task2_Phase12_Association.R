# ==========================================
# Task 2 - Phase 12: Metabolite-Diabetes Association
# ==========================================

cat("=== Phase 12: Metabolite-Diabetes Phenotype Association ===\n\n")

if (!require("dplyr", quietly=TRUE)) install.packages("dplyr")
if (!require("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!require("tidyr", quietly=TRUE)) install.packages("tidyr")
if (!require("ggrepel", quietly=TRUE)) install.packages("ggrepel")

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)

# 1. Load Data
cat("Loading datasets...\n")
qc_data <- read.csv("../intermediate/Phase11_Metabolites_QC_Final.csv", stringsAsFactors = FALSE)

# Load covariates from Task 1
# Sex
sex_df <- read.table("../intermediate/pheno_sex.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(sex_df) <- c("FID", "main_id", "Sex")
# PC1, PC2
pc1_df <- read.table("../intermediate/pheno_PC1.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(pc1_df) <- c("FID", "main_id", "PC1")
pc2_df <- read.table("../intermediate/pheno_PC2.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(pc2_df) <- c("FID", "main_id", "PC2")
# PC3-10
covar_df <- read.table("../intermediate/covar.txt", header=TRUE, stringsAsFactors=FALSE)
colnames(covar_df)[1:2] <- c("FID", "main_id")

# Merge covariates
covariates <- sex_df %>%
  left_join(pc1_df, by=c("FID", "main_id")) %>%
  left_join(pc2_df, by=c("FID", "main_id")) %>%
  left_join(covar_df, by=c("FID", "main_id"))

# Merge covariates into qc_data
full_data <- qc_data %>%
  left_join(covariates, by="main_id")

# Convert Diabetes to factor for plotting
full_data$Diabetes_Status <- factor(full_data$Diabetes, levels=c(0, 1), labels=c("Non-Diabetic", "Diabetic"))

# 2. Exploratory Bar Plot
cat("Generating exploratory bar plot...\n")

metabolite_cols <- names(qc_data)[!(names(qc_data) %in% c("main_id", "mapped_id", "Diabetes"))]

summary_stats <- full_data %>%
  select(Diabetes_Status, all_of(metabolite_cols)) %>%
  pivot_longer(cols = -Diabetes_Status, names_to = "Metabolite", values_to = "Level") %>%
  group_by(Metabolite, Diabetes_Status) %>%
  summarise(
    Mean = mean(Level, na.rm = TRUE),
    SE = sd(Level, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Since 136 metabolites is too many for one plot, we will plot the top 40 by absolute difference in means
mean_diffs <- summary_stats %>%
  select(Metabolite, Diabetes_Status, Mean) %>%
  pivot_wider(names_from = Diabetes_Status, values_from = Mean) %>%
  mutate(AbsDiff = abs(Diabetic - `Non-Diabetic`)) %>%
  arrange(desc(AbsDiff))

top_metabolites <- head(mean_diffs$Metabolite, 40)
plot_data <- summary_stats %>% filter(Metabolite %in% top_metabolites)

pdf("../outputs/Phase12_Association/Phase12_Exploratory_BarPlot.pdf", width=12, height=10)
print(
  ggplot(plot_data, aes(x = reorder(Metabolite, Mean), y = Mean, fill = Diabetes_Status)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(0.9), width = 0.25) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "Metabolite Levels by Diabetes Status (Top 40 by Mean Difference)",
      x = "Metabolite",
      y = "Mean Z-Score Level (± SE)",
      fill = "Status"
    ) +
    scale_fill_manual(values = c("Non-Diabetic" = "steelblue", "Diabetic" = "darkred"))
)
invisible(dev.off())

# 3. Formal Association Test (Logistic Regression)
cat("Running formal association tests (logistic regression)...\n")

results_list <- list()

# Prepare model data, removing rows with missing covariates if any
model_data <- full_data %>% drop_na(Diabetes, Sex, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

for (metab in metabolite_cols) {
  # Formula: Diabetes ~ Metabolite + Sex + PC1 ... PC10
  f <- as.formula(paste("Diabetes ~", metab, "+ Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
  
  # Fit logistic regression
  fit <- glm(f, data = model_data, family = binomial)
  
  # Extract summary
  summ <- summary(fit)
  
  # Check if metabolite is in the summary (sometimes dropped if perfectly collinear, though rare here)
  if (metab %in% rownames(summ$coefficients)) {
    beta <- summ$coefficients[metab, "Estimate"]
    se <- summ$coefficients[metab, "Std. Error"]
    p_val <- summ$coefficients[metab, "Pr(>|z|)"]
    
    results_list[[metab]] <- data.frame(
      Metabolite = metab,
      Beta = beta,
      SE = se,
      P_value = p_val,
      stringsAsFactors = FALSE
    )
  }
}

assoc_results <- do.call(rbind, results_list)

# 4. Bonferroni Correction
cat("Applying Bonferroni correction...\n")
n_tests <- nrow(assoc_results)
alpha <- 0.05
bonferroni_thresh <- alpha / n_tests

assoc_results <- assoc_results %>%
  mutate(
    Bonferroni_Threshold = bonferroni_thresh,
    Significant = ifelse(P_value < bonferroni_thresh, "Y", "N")
  ) %>%
  arrange(P_value)

write.csv(assoc_results, "../outputs/Phase12_Association/Phase12_Metabolite_Association_Results.csv", row.names = FALSE)

num_sig <- sum(assoc_results$Significant == "Y")
cat(sprintf("Found %d significant metabolites (p < %.2e).\n", num_sig, bonferroni_thresh))

# Save the list of significant metabolites for Phase 13
sig_metabs <- assoc_results$Metabolite[assoc_results$Significant == "Y"]
write.table(sig_metabs, "../intermediate/Phase12_Significant_Metabolites.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 5. Volcano Plot
cat("Generating Volcano Plot...\n")
assoc_results <- assoc_results %>%
  mutate(
    log10P = -log10(P_value),
    Plot_Label = ifelse(Significant == "Y", Metabolite, "")
  )

pdf("../outputs/Phase12_Association/Phase12_Volcano_Plot.pdf", width=10, height=8)
print(
  ggplot(assoc_results, aes(x = Beta, y = log10P, color = Significant)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = -log10(bonferroni_thresh), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_text_repel(aes(label = Plot_Label), size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("Y" = "darkred", "N" = "darkgray")) +
    theme_minimal() +
    labs(
      title = "Volcano Plot: Metabolite-Diabetes Association",
      subtitle = sprintf("Logistic regression adjusted for Sex and PC1-10\nRed dashed line = Bonferroni threshold (p < %.2e)", bonferroni_thresh),
      x = "Log Odds Ratio (Beta)",
      y = "-log10(P-value)"
    )
)
invisible(dev.off())

cat("\nPhase 12 complete. Outputs saved to data/:\n")
cat("  - Phase12_Exploratory_BarPlot.pdf\n")
cat("  - Phase12_Metabolite_Association_Results.csv\n")
cat("  - Phase12_Volcano_Plot.pdf\n")
cat("  - intermediate_files/Phase12_Significant_Metabolites.txt\n")
