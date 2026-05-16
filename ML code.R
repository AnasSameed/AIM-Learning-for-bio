# ============================================
# GSE308050 Analysis (Updated)
# Description: Methylation Analysis & Random Forest Classifier
# ============================================

# Load Required Libraries
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(pheatmap)
library(dplyr)

set.seed(42)

# 1. SETUP PATHS
data_path <- "C:/Users/User/Desktop/Dr. Aniruddha/Apr-May learnings/GSE308050/GSE308050_All_merged_DMRs_sig_methy_group_clinics_infor.csv.gz"
plot_path <- "C:/Users/User/Desktop/Dr. Aniruddha/Apr-May learnings/GSE308050/plots/"
res_path  <- "C:/Users/User/Desktop/Dr. Aniruddha/Apr-May learnings/GSE308050/results/"

dir.create(plot_path, showWarnings = FALSE)
dir.create(res_path,  showWarnings = FALSE)

# ============================================
# STEP 1: Load Data
# ============================================
cat("Loading data...\n")
df_raw <- read.csv(gzfile(data_path), header = TRUE, row.names = 1, check.names = FALSE)
cat("Raw dimensions:", dim(df_raw), "\n")

# ============================================
# STEP 2: Separate Clinical and Methylation Data
# ============================================
clinical_vars <- c("surg_pgleasnt", "surg_pstagt", "clinicalt",
                   "recurrence", "Death", "Bad_outcome",
                   "PSA_recurrence", "surg_lnmets", "surg_age", "risky")

clinical             <- as.data.frame(t(df_raw[clinical_vars, ]))
clinical$Bad_outcome <- factor(ifelse(clinical$Bad_outcome == 1, "Bad", "Good"), levels = c("Good", "Bad"))
clinical$risky       <- factor(clinical$risky, levels = c("Low","Middle","High"))
clinical$surg_age    <- as.numeric(clinical$surg_age)

dmr_rows   <- setdiff(rownames(df_raw), clinical_vars)
methyl_raw <- df_raw[dmr_rows, ]
methyl     <- as.data.frame(t(apply(methyl_raw, 1, as.numeric)))
colnames(methyl) <- colnames(df_raw)
methyl     <- as.data.frame(t(methyl))

cat("Class distribution:\n")
print(table(clinical$Bad_outcome))

# ============================================
# STEP 3: Feature Selection (Spearman Correlation)
# ============================================
cat("\nSelecting top 200 features based on correlation...\n")
outcome_numeric <- as.numeric(clinical$Bad_outcome == "Bad")
cors <- apply(methyl, 2, function(x) cor(x, outcome_numeric, method = "spearman", use = "complete.obs"))
top_dmrs <- names(sort(abs(cors), decreasing = TRUE)[1:200])

ml_data       <- as.data.frame(methyl[, top_dmrs])
ml_data$Label <- clinical$Bad_outcome

# ============================================
# STEP 4: Cross-Validated Random Forest
# ============================================
cat("\nRunning 5-fold repeated cross-validation (5x5)...\n")

ctrl <- trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats          = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Calculate class weights for imbalance
class_weights <- c(
  Bad  = nrow(ml_data) / (2 * sum(ml_data$Label == "Bad")),
  Good = nrow(ml_data) / (2 * sum(ml_data$Label == "Good"))
)

rf_model <- train(
  Label ~ .,
  data       = ml_data,
  method     = "rf",
  metric     = "ROC",
  trControl  = ctrl,
  importance = TRUE,
  classwt    = class_weights,
  tuneGrid   = expand.grid(mtry = c(5, 10, 15, 20))
)

cat("\nBest mtry:", rf_model$bestTune$mtry, "\n")
cat("Final CV AUC:", round(max(rf_model$results$ROC), 3), "\n")

# ============================================
# STEP 5: ROC Curve Generation
# ============================================
cv_preds <- rf_model$pred
roc_obj  <- roc(cv_preds$obs, cv_preds$Bad, levels = c("Good","Bad"))

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "#e74c3c", size = 1.2) +
  geom_abline(linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.6, y = 0.2,
           label = paste0("AUC = ", round(auc(roc_obj), 3)),
           size = 5, color = "#2c3e50", fontface = "bold") +
  labs(title    = "ROC Curve — Random Forest (Cross-Validated)",
       subtitle = "Predicting Bad Outcome from DMR Methylation",
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal(base_size = 13)

ggsave(paste0(plot_path, "ROC_curve_CV.png"), p_roc, width = 7, height = 6, dpi = 300)

# ============================================
# STEP 6: Robust Feature Importance (The Fix)
# ============================================
cat("\nGenerating Feature Importance Plot...\n")

imp_raw <- varImp(rf_model)$importance
imp_df  <- as.data.frame(imp_raw)
imp_df$DMR <- rownames(imp_df)

# Check for column types and assign Score
if ("Overall" %in% colnames(imp_df)) {
  imp_df$Score <- imp_df$Overall
} else if ("Bad" %in% colnames(imp_df)) {
  imp_df$Score <- imp_df$Bad
} else {
  imp_df$Score <- rowMeans(imp_df[, sapply(imp_df, is.numeric)])
}

# Sort and take top 20
imp_df <- imp_df[order(imp_df$Score, decreasing = TRUE), ]
top20  <- head(imp_df, 20)
top20$DMR <- factor(top20$DMR, levels = rev(top20$DMR)) # For correct ggplot ordering

p_imp <- ggplot(top20, aes(x = DMR, y = Score)) +
  geom_bar(stat = "identity", fill = "#2c3e50", alpha = 0.85) +
  coord_flip() +
  labs(title    = "Top 20 Predictive DMRs",
       subtitle = "Random Forest Importance Score",
       x = "DMR Region", y = "Importance Score") +
  theme_minimal(base_size = 11)

ggsave(paste0(plot_path, "feature_importance_CV.png"), p_imp, width = 10, height = 7, dpi = 300)
write.csv(imp_df, paste0(res_path, "DMR_importance_ranked.csv"), row.names = FALSE)

cat("\nAnalysis Complete. Files saved in 'plots/' and 'results/' folders.\n")