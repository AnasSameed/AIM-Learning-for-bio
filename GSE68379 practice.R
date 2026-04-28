# ============================================
# Breast Cancer DNA Methylation Analysis
# Dataset: GSE68379 (GEO)
# Learning....
# ============================================

# STEP 0: Install packages (run once only)
# install.packages("BiocManager")
# BiocManager::install(c("GEOquery", "limma"))
# install.packages("ggplot2")
# install.packages("pheatmap")

library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)

# ============================================
# STEP 1: Download data from GEO
# ============================================
cat("Downloading GSE68379...\n")
gse <- getGEO("GSE68379", GSEMatrix = TRUE, getGPL = FALSE)
gse_data <- gse[[1]]

# Extract beta values and phenotype
beta_matrix <- exprs(gse_data)
pheno_data  <- pData(gse_data)

cat("Dataset dimensions:", dim(beta_matrix), "\n")
cat("Samples:", ncol(beta_matrix), "\n")
cat("CpG sites:", nrow(beta_matrix), "\n")

# ============================================
# STEP 2: Check sample groups
# ============================================
cat("\n--- Sample Info ---\n")
print(pheno_data[, c("title", "source_name_ch1")])

# ============================================
# STEP 3: Filter CpGs
# ============================================
cat("\nFiltering CpGs...\n")

# Remove rows with any NA
beta_matrix <- beta_matrix[complete.cases(beta_matrix), ]

# Keep top 25% most variable CpGs
variance <- apply(beta_matrix, 1, var)
top_var  <- beta_matrix[variance > quantile(variance, 0.75), ]

cat("CpGs after filtering:", nrow(top_var), "\n")

# ============================================
# STEP 4: PCA Plot
# ============================================
cat("\nGenerating PCA plot...\n")
dir.create("plots", showWarnings = FALSE)

pca_result <- prcomp(t(top_var), scale. = TRUE)
pca_df <- data.frame(
  PC1    = pca_result$x[, 1],
  PC2    = pca_result$x[, 2],
  Sample = rownames(pca_result$x)
)

ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(vjust = -0.5, size = 2.5) +
  labs(
    title = "PCA — Breast Cancer Methylation (GSE68379)",
    x = "PC1", y = "PC2"
  ) +
  theme_minimal()

ggsave("plots/PCA_plot.png", width = 8, height = 6)
cat("PCA plot saved.\n")

# ============================================
# STEP 5: Define groups (tumor vs normal)
# ============================================
# After running Step 2, look at source_name_ch1
# Then define groups like below (adjust labels to match your data)

# Example — update after checking pheno_data output:
# group <- factor(ifelse(grepl("tumor", pheno_data$source_name_ch1,
#                              ignore.case = TRUE), "Tumor", "Normal"))

# ============================================
# STEP 6: Differential Methylation (limma)
# ============================================
# Run this AFTER defining group above

# design  <- model.matrix(~group)
# fit     <- lmFit(top_var, design)
# fit     <- eBayes(fit)
# results <- topTable(fit, coef = 2, number = Inf)
# results <- results[order(results$adj.P.Val), ]

# dir.create("results", showWarnings = FALSE)
# write.csv(results, "results/differential_methylation.csv", row.names = TRUE)
# cat("Differential methylation results saved.\n")

# ============================================
# STEP 7: Heatmap (top 50 DMPs)
# ============================================
# Run after Step 6

# top50 <- head(rownames(results), 50)
# pheatmap(
#   top_var[top50, ],
#   show_rownames    = FALSE,
#   show_colnames    = TRUE,
#   scale            = "row",
#   main             = "Top 50 Differentially Methylated CpGs",
#   filename         = "plots/heatmap_top50.png"
# )
# cat("Heatmap saved.\n")

cat("\nSteps 1-4 complete. Check plots/ folder for PCA.\n")
cat("Paste your pheno_data output to proceed with Step 5.\n")