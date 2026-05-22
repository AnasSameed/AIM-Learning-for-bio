
# Prostate Cancer Methylation Analysis
# Dataset: GSE308050

#Load lib

library(ggplot2)
library(dplyr)

# STEP 1: Load data

df_raw <- read.csv(
  gzfile("C:/Users/User/Desktop/Dr. Aniruddha/Apr-May learnings/GSE308050/GSE308050_All_merged_DMRs_sig_methy_group_clinics_infor.csv.gz"),
  header = TRUE, row.names = 1, check.names = FALSE
)

cat("Raw dimensions:", dim(df_raw), "\n")
cat("First few row names:\n")
print(rownames(df_raw)[1:15])

# STEP 2: Separate clinical and methylation data

clinical_vars <- c("surg_pgleasnt", "surg_pstagt", "clinicalt",
                   "recurrence", "Death", "Bad_outcome",
                   "PSA_recurrence", "surg_lnmets", "surg_age", "risky")

# Clinical matrix (10 variables x 120 patients) → transpose to patients x variables
clinical <- as.data.frame(t(df_raw[clinical_vars, ]))
clinical$Bad_outcome  <- as.factor(clinical$Bad_outcome)
clinical$risky        <- factor(clinical$risky, levels = c("Low","Middle","High"))
clinical$surg_age     <- as.numeric(clinical$surg_age)
clinical$surg_pgleasnt <- as.numeric(clinical$surg_pgleasnt)

# Methylation matrix → transpose to patients x DMRs
dmr_rows   <- setdiff(rownames(df_raw), clinical_vars)
methyl_raw <- df_raw[dmr_rows, ]
methyl     <- as.data.frame(t(apply(methyl_raw, 1, as.numeric)))
colnames(methyl) <- colnames(df_raw)
methyl     <- as.data.frame(t(methyl))  # now patients x DMRs

cat("Clinical data:", dim(clinical), "\n")
cat("Methylation data:", dim(methyl), "\n")

# STEP 3: Quick summary

cat("\n--- Risk group distribution ---\n")
print(table(clinical$risky))

cat("\n--- Bad outcome distribution ---\n")
print(table(clinical$Bad_outcome))

cat("\n--- Age summary ---\n")
print(summary(clinical$surg_age))

# STEP 4: PCA on methylation data

dir.create("C:/Users/User/Desktop/Dr. Aniruddha/Apr-May learnings/GSE308050/plots", 
           showWarnings = FALSE)
dir.create("C:/Users/User/Desktop/Dr. Aniruddha/Apr-May learnings/GSE308050/results", 
           showWarnings = FALSE)

plot_path <- "C:/Users/User/Desktop/Dr. Aniruddha/Apr-May learnings/GSE308050/plots/"

# Remove zero variance columns before PCA
var_filter <- apply(methyl, 2, var, na.rm = TRUE)
methyl_filt <- methyl[, var_filter > 0]
cat("DMRs after variance filter:", ncol(methyl_filt), "\n")

# Run PCA
pca_result <- prcomp(methyl_filt, scale. = TRUE)
pca_df <- data.frame(
  PC1   = pca_result$x[, 1],
  PC2   = pca_result$x[, 2],
  PC3   = pca_result$x[, 3],
  Risk  = clinical$risky,
  Outcome = clinical$Bad_outcome,
  Age   = clinical$surg_age
)

# PCA colored by Risk group
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Risk)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Low" = "#2ecc71", 
                                "Middle" = "#f39c12", 
                                "High" = "#e74c3c")) +
  labs(title = "PCA of DNA Methylation — Colored by Risk Group",
       subtitle = "GSE308050 | 120 Prostate Cancer Patients",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "% variance)")) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(paste0(plot_path, "PCA_by_risk.png"), p1, width = 8, height = 6, dpi = 300)

# PCA colored by Bad outcome
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Outcome)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#3498db", "1" = "#e74c3c"),
                     labels = c("0" = "No Bad Outcome", "1" = "Bad Outcome")) +
  labs(title = "PCA of DNA Methylation — Colored by Outcome",
       subtitle = "GSE308050 | 120 Prostate Cancer Patients",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "% variance)")) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(paste0(plot_path, "PCA_by_outcome.png"), p2, width = 8, height = 6, dpi = 300)

cat("PCA plots saved.\n")

# STEP 5: Variance explained plot

var_explained <- summary(pca_result)$importance[2, 1:20] * 100
ve_df <- data.frame(
  PC       = paste0("PC", 1:20),
  Variance = var_explained
)
ve_df$PC <- factor(ve_df$PC, levels = ve_df$PC)

p3 <- ggplot(ve_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#2c3e50") +
  geom_line(aes(group = 1), color = "#e74c3c", size = 1) +
  geom_point(color = "#e74c3c", size = 2) +
  labs(title = "Variance Explained by Top 20 Principal Components",
       x = "Principal Component", y = "% Variance Explained") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(paste0(plot_path, "PCA_variance_explained.png"), p3, width = 9, height = 5, dpi = 300)
cat("Variance plot saved.\n")

# STEP 6: Top variable DMRs heatmap

library(pheatmap)

# Top 50 most variable DMRs
top_var_idx  <- order(var_filter, decreasing = TRUE)[1:50]
top50_methyl <- methyl[, top_var_idx]

# Annotation for heatmap
annot_df <- data.frame(
  Risk    = clinical$risky,
  Outcome = ifelse(clinical$Bad_outcome == 1, "Bad", "Good"),
  row.names = rownames(clinical)
)

annot_colors <- list(
  Risk    = c(Low = "#2ecc71", Middle = "#f39c12", High = "#e74c3c"),
  Outcome = c(Bad = "#e74c3c", Good = "#3498db")
)

pheatmap(
  t(top50_methyl),
  annotation_col  = annot_df,
  annotation_colors = annot_colors,
  show_colnames   = FALSE,
  show_rownames   = FALSE,
  scale           = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main            = "Top 50 Variable DMRs — Methylation Heatmap",
  filename        = paste0(plot_path, "heatmap_top50.png"),
  width = 10, height = 8
)
cat("Heatmap saved.\n")

# ============================================
# STEP 7: Risk group methylation distribution
# ============================================
# Average methylation per patient
avg_methyl <- rowMeans(methyl, na.rm = TRUE)

dist_df <- data.frame(
  AvgMethylation = avg_methyl,
  Risk = clinical$risky,
  Outcome = clinical$Bad_outcome
)

p4 <- ggplot(dist_df, aes(x = Risk, y = AvgMethylation, fill = Risk)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = c("Low" = "#2ecc71", 
                               "Middle" = "#f39c12", 
                               "High" = "#e74c3c")) +
  labs(title = "Average DMR Methylation by Risk Group",
       x = "Risk Group", y = "Average Methylation (beta)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")

ggsave(paste0(plot_path, "methylation_by_risk.png"), p4, width = 7, height = 6, dpi = 300)
cat("Distribution plot saved.\n")

cat("\nAll Step 4-7 plots saved to plots/ folder.\n")
cat("Open the plots folder and share screenshots of the plots!\n")

# ============================================
# STEP 8: Random Forest ML Classifier
# ============================================
install.packages("randomForest")
install.packages("caret")
install.packages("pROC")

#library(randomForest)
#library(caret)
#library(pROC)

set.seed(42)

# Prepare data — use top 500 most variable DMRs for speed
top500_idx    <- order(var_filter, decreasing = TRUE)[1:500]
ml_data       <- methyl[, top500_idx]
ml_data$Label <- clinical$Bad_outcome

# Train/test split — 80/20
train_idx <- createDataPartition(ml_data$Label, p = 0.8, list = FALSE)
train_df  <- ml_data[train_idx, ]
test_df   <- ml_data[-train_idx, ]

cat("Training samples:", nrow(train_df), "\n")
cat("Testing samples:", nrow(test_df), "\n")

# Train Random Forest
cat("Training Random Forest...\n")
rf_model <- randomForest(
  Label ~ .,
  data       = train_df,
  ntree      = 500,
  importance = TRUE,
  mtry       = 22  # sqrt(500)
)

print(rf_model)

# Predict on test set
predictions  <- predict(rf_model, test_df, type = "prob")[, 2]
pred_class   <- predict(rf_model, test_df)
actual       <- test_df$Label

# Confusion matrix
cat("\n--- Confusion Matrix ---\n")
print(confusionMatrix(pred_class, actual))

# ROC curve
roc_obj <- roc(as.numeric(actual), predictions)
auc_val <- auc(roc_obj)
cat("\nAUC:", round(auc_val, 3), "\n")

# Plot ROC
roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "#e74c3c", size = 1.2) +
  geom_abline(linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.65, y = 0.15,
           label = paste0("AUC = ", round(auc_val, 3)),
           size = 5, color = "#2c3e50", fontface = "bold") +
  labs(title = "ROC Curve — Random Forest Classifier",
       subtitle = "Predicting Bad Outcome from DMR Methylation",
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(paste0(plot_path, "ROC_curve.png"), p_roc, width = 7, height = 6, dpi = 300)

# Feature importance — top 20 DMRs
imp_df <- as.data.frame(importance(rf_model))
imp_df$DMR <- rownames(imp_df)
imp_df <- imp_df[order(imp_df$MeanDecreaseGini, decreasing = TRUE), ]
top20_imp <- head(imp_df, 20)

p_imp <- ggplot(top20_imp, aes(x = reorder(DMR, MeanDecreaseGini),
                               y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#2c3e50") +
  coord_flip() +
  labs(title = "Top 20 Most Important DMRs",
       subtitle = "Random Forest Feature Importance (Mean Decrease Gini)",
       x = "DMR Region", y = "Importance Score") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(paste0(plot_path, "feature_importance.png"), p_imp,
       width = 9, height = 7, dpi = 300)

# Save top 20 important DMRs
write.csv(top20_imp, paste0(
  "C:/Users/User/Desktop/Dr. Aniruddha/Apr-May learnings/GSE308050/results/top20_important_DMRs.csv"),
  row.names = FALSE)

cat("\nML analysis complete. Check plots/ for ROC and importance plots.\n")