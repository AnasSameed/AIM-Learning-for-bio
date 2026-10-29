# TCGA-PRAD tumor vs normal: RNA-seq differential expression pipeline (R)
# Mirrors my TCGA-PRAD analysis report. TEMPLATE: re-run and check versions before reuse.
# install: BiocManager::install(c("TCGAbiolinks","DESeq2","clusterProfiler","org.Hs.eg.db"))
#          install.packages(c("ggplot2","pheatmap","survminer","survival","ggrepel"))
suppressPackageStartupMessages({
  library(TCGAbiolinks); library(SummarizedExperiment); library(DESeq2)
  library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)
  library(pheatmap); library(ggrepel)
})

# 1. Download STAR counts from GDC
q <- GDCquery(project = "TCGA-PRAD", data.category = "Transcriptome Profiling",
              data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(q); se <- GDCprepare(q)

# 2. Subset 50 tumor + 10 normal
set.seed(1)
tum <- sample(which(se$sample_type == "Primary Tumor"), 50)
nor <- sample(which(se$sample_type == "Solid Tissue Normal"), 10)
se  <- se[, c(tum, nor)]

# 3. DESeq2 (filter low counts, BH-adjusted p)
dds <- DESeqDataSetFromMatrix(assay(se, "unstranded"), colData(se), design = ~ sample_type)
dds$sample_type <- relevel(factor(dds$sample_type), ref = "Solid Tissue Normal")
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
dds <- DESeq(dds)
res <- as.data.frame(results(dds, alpha = 0.05)); res$gene <- rownames(res)
sig <- subset(res, abs(log2FoldChange) > 1 & padj < 0.05)
cat("Significant DEGs:", nrow(sig), "\n")

# 4. Volcano plot
res$status <- ifelse(res$padj < 0.05 & res$log2FoldChange >  1, "Up",
              ifelse(res$padj < 0.05 & res$log2FoldChange < -1, "Down", "NS"))
ggplot(res, aes(log2FoldChange, -log10(padj), colour = status)) + geom_point(size = 0.7) +
  scale_colour_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) + geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_bw() + ggtitle("Volcano: TCGA-PRAD tumor vs normal")
ggsave("volcano.png", width = 7, height = 5, dpi = 300)

# 5. PCA and heatmap of top 50 genes
vsd <- vst(dds, blind = TRUE)
ggsave("pca.png", plotPCA(vsd, intgroup = "sample_type") + theme_bw(), width = 6, height = 5)
top <- head(rownames(res[order(res$padj), ]), 50)
pheatmap(assay(vsd)[top, ], scale = "row", show_rownames = FALSE,
         annotation_col = as.data.frame(colData(dds)["sample_type"]), filename = "heatmap_top50.png")

# 6. GO enrichment (up and down separately)
to_entrez <- function(ens) bitr(sub("\\..*", "", ens), "ENSEMBL", "ENTREZID", org.Hs.eg.db)$ENTREZID
up   <- to_entrez(rownames(subset(sig, log2FoldChange > 0)))
down <- to_entrez(rownames(subset(sig, log2FoldChange < 0)))
ggsave("go_up.png",   dotplot(enrichGO(up,   org.Hs.eg.db, ont = "BP", readable = TRUE)), width = 7, height = 6)
ggsave("go_down.png", dotplot(enrichGO(down, org.Hs.eg.db, ont = "BP", readable = TRUE)), width = 7, height = 6)

# 7. Export a VST matrix + labels for the Python ML pipeline
write.csv(t(assay(vsd)), "vst_matrix.csv")
write.csv(data.frame(sample = colnames(dds), group = dds$sample_type), "labels.csv", row.names = FALSE)
