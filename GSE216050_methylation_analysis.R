############################################################
#  COMPREHENSIVE RRBS DNA METHYLATION ANALYSIS
#  Dataset : GSE216050  (Prostate Cancer — Decitabine Treatment)
#  3 Controls (PM154, PM12, PM21) vs 3 Decitabine-treated
#  Supervisor: Dr Aniruddha Chatterjee | University of Otago
#
#  ANALYSES COVERED:
#    1.  Data loading & quality control
#    2.  Global methylation comparison (boxplot + violin)
#    3.  CpG coverage diagnostics
#    4.  PCA of methylation profiles
#    5.  Correlation heatmap between samples
#    6.  Differential methylation (t-test + FDR correction)
#    7.  Volcano plot (publication-quality)
#    8.  DMR identification via sliding window
#    9.  Genomic annotation of DMRs (promoter/exon/intron/intergenic)
#   10.  CpG island context annotation
#   11.  Heatmap of top differentially methylated CpGs
#   12.  GO enrichment analysis (Biological Process)
#   13.  KEGG pathway enrichment
#   14.  Methylation density plots per chromosome
#   15.  B7-H3 (CD276) locus-specific methylation plot
#   16.  Summary report table (CSV)
############################################################

# ─────────────────────────────────────────────────────────
#  0. INSTALL / LOAD PACKAGES
# ─────────────────────────────────────────────────────────

# Run the install block once; comment out afterwards
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c(
#   "GenomicRanges","IRanges","GenomicFeatures",
#   "TxDb.Hsapiens.UCSC.hg38.knownGene",
#   "org.Hs.eg.db","clusterProfiler",
#   "AnnotationHub","ChIPseeker","annotatr"
# ))
# install.packages(c("data.table","ggplot2","dplyr","pheatmap",
#                    "ggrepel","RColorBrewer","viridis",
#                    "patchwork","scales","tidyr"))
#install.packages(c(
#  "viridis", "patchwork", "scales",
#  "ggrepel", "RColorBrewer", "data.table",
#  "ggplot2", "dplyr", "tidyr", "pheatmap"
#))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(patchwork)
  library(scales)
  library(GenomicRanges)
  library(IRanges)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

# ─────────────────────────────────────────────────────────
#  GLOBAL SETTINGS
# ─────────────────────────────────────────────────────────
setwd("C:/Users/User/Downloads/demo/Methylation dataset/GSE216050_RAW")
set.seed(42)

# Output directory — creates a dated subfolder
out_dir <- file.path("results", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Thresholds — adjust here to explore sensitivity
MIN_COVERAGE   <- 10      # minimum read depth per CpG
DIFF_THRESHOLD <- 20      # absolute methylation difference (%)
PVAL_CUTOFF    <- 0.05    # raw p-value cutoff
FDR_CUTOFF     <- 0.05    # BH-adjusted FDR cutoff
DMR_WIN_SIZE   <- 1000    # sliding window size for DMR calling (bp)
DMR_MIN_CPGS   <- 3       # minimum CpGs per DMR

# Sample metadata — edit paths if files are in a subfolder
sample_meta <- data.frame(
  sample_id  = c("C1","C2","C3","T1","T2","T3"),
  file       = c(
    "GSM6657904_data1.txt",
    "GSM6657906_data3.txt",
    "GSM6657908_data7.txt",
    "GSM6657905_data2.txt",
    "GSM6657907_data5.txt",
    "GSM6657909_data8.txt"
  ),
  group      = c(rep("Control",3), rep("Decitabine",3)),
  patient    = c("PM154","PM12","PM21","PM154","PM12","PM21"),
  stringsAsFactors = FALSE
)

# Colour palette
group_colours <- c(Control = "#2166AC", Decitabine = "#D6604D")
sample_colours <- c(
  C1="#4393C3", C2="#2166AC", C3="#053061",
  T1="#F4A582", T2="#D6604D", T3="#67001F"
)

cat("============================================================\n")
cat(" GSE216050 — RRBS Methylation Analysis Pipeline\n")
cat(" Output directory:", out_dir, "\n")
cat("============================================================\n\n")

# ─────────────────────────────────────────────────────────
#  STEP 1 — READ & QC INDIVIDUAL SAMPLE FILES
# ─────────────────────────────────────────────────────────

cat("[1/16] Reading and quality-controlling sample files...\n")

read_rrbs <- function(file, sample_name, min_cov = MIN_COVERAGE) {
  
  if (!file.exists(file)) {
    stop("File not found: ", file,
         "\n  Please check your working directory: ", getwd())
  }
  
  dt <- fread(file, header = FALSE, showProgress = FALSE)
  
  # Handle both 5-column and 6-column bismark coverage2cytosine output
  if (ncol(dt) == 5) {
    setnames(dt, c("chr","start","end","meth_pct","meth_reads"))
    dt[, unmeth_reads := NA_integer_]
  } else {
    setnames(dt, c("chr","start","end","meth_pct","meth_reads","unmeth_reads"))
  }
  
  # Compute coverage
  dt[, coverage := meth_reads + unmeth_reads]
  
  # Basic QC metrics BEFORE filtering
  qc <- list(
    sample        = sample_name,
    total_cpgs    = nrow(dt),
    cpgs_cov10    = sum(dt$coverage >= min_cov, na.rm = TRUE),
    median_cov    = median(dt$coverage, na.rm = TRUE),
    mean_meth     = mean(dt$meth_pct, na.rm = TRUE),
    pct_hyper     = mean(dt$meth_pct >= 80, na.rm = TRUE) * 100,
    pct_hypo      = mean(dt$meth_pct <= 20, na.rm = TRUE) * 100
  )
  
  # Filter by minimum coverage
  dt <- dt[coverage >= min_cov]
  
  # Keep only autosomes + chrX for cleaner analysis
  keep_chr <- c(paste0("chr", 1:22), "chrX")
  dt <- dt[chr %in% keep_chr]
  
  # Normalise methylation to 0–1 (beta value) alongside %
  dt[, beta := meth_pct / 100]
  
  # Return tidy per-CpG table
  dt[, .(chr, start, end,
         meth_pct = round(meth_pct, 4),
         beta     = round(beta, 4),
         coverage)]
}

# Read all samples
sample_list <- mapply(
  read_rrbs,
  file        = sample_meta$file,
  sample_name = sample_meta$sample_id,
  SIMPLIFY    = FALSE
)
names(sample_list) <- sample_meta$sample_id

# Collect QC summary
qc_list <- lapply(seq_along(sample_list), function(i) {
  s <- sample_list[[i]]
  data.frame(
    Sample       = names(sample_list)[i],
    Group        = sample_meta$group[i],
    Total_CpGs   = nrow(s),
    Median_Cov   = median(s$coverage),
    Mean_Meth_pct = round(mean(s$meth_pct), 2),
    stringsAsFactors = FALSE
  )
})
qc_df <- do.call(rbind, qc_list)
write.csv(qc_df, file.path(out_dir, "00_QC_summary.csv"), row.names = FALSE)
print(qc_df)

# ─────────────────────────────────────────────────────────
#  STEP 2 — MERGE ON COMMON CpG SITES
# ─────────────────────────────────────────────────────────

cat("\n[2/16] Merging samples on common CpG sites...\n")

# Rename columns before merging
for (s in names(sample_list)) {
  setnames(sample_list[[s]], "meth_pct", s)
  setnames(sample_list[[s]], "beta", paste0(s, "_beta"))
  sample_list[[s]][, coverage := NULL]   # drop per-sample coverage to save RAM
}

merged <- Reduce(
  function(x, y) merge(x, y, by = c("chr","start","end")),
  sample_list
)

cat("   Common CpG sites across all 6 samples:", format(nrow(merged), big.mark = ","), "\n")

# ─────────────────────────────────────────────────────────
#  STEP 3 — COVERAGE DISTRIBUTION PLOT
# ─────────────────────────────────────────────────────────

cat("[3/16] Plotting per-sample coverage distributions...\n")

# Re-read with coverage for plotting
cov_plot_data <- lapply(seq_len(nrow(sample_meta)), function(i) {
  dt <- fread(sample_meta$file[i], header = FALSE, showProgress = FALSE)
  if (ncol(dt) == 6)
    setnames(dt, c("chr","start","end","meth_pct","meth_reads","unmeth_reads"))
  else
    setnames(dt, c("chr","start","end","meth_pct","meth_reads","unmeth_reads"))
  dt[, cov := meth_reads + unmeth_reads]
  dt[cov >= 1 & cov <= 100,
     .(coverage = cov, sample = sample_meta$sample_id[i],
       group = sample_meta$group[i])]
})
cov_df <- rbindlist(cov_plot_data)

p_cov <- ggplot(cov_df, aes(x = coverage, fill = group)) +
  geom_histogram(binwidth = 2, alpha = 0.7, position = "identity") +
  facet_wrap(~sample, nrow = 2) +
  scale_fill_manual(values = group_colours) +
  scale_x_continuous(limits = c(0, 60)) +
  theme_classic(base_size = 12) +
  labs(title = "CpG Coverage Distribution per Sample",
       subtitle = "GSE216050 | RRBS | Prostate Cancer",
       x = "Read Depth", y = "Number of CpGs") +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "01_Coverage_distribution.png"),
       p_cov, width = 12, height = 6, dpi = 300)

# ─────────────────────────────────────────────────────────
#  STEP 4 — GLOBAL METHYLATION COMPARISON
# ─────────────────────────────────────────────────────────

cat("[4/16] Global methylation comparison plots...\n")

long_meth <- merged %>%
  select(chr, start, end, C1, C2, C3, T1, T2, T3) %>%
  pivot_longer(cols = C1:T3, names_to = "sample", values_to = "methylation") %>%
  left_join(sample_meta[, c("sample_id","group","patient")],
            by = c("sample" = "sample_id"))

# Violin + boxplot overlay
p_violin <- ggplot(long_meth,
                   aes(x = sample, y = methylation, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.9) +
  scale_fill_manual(values = group_colours) +
  theme_classic(base_size = 13) +
  labs(title = "Global CpG Methylation — Control vs Decitabine",
       subtitle = "Common CpGs (coverage ≥ 10)",
       x = NULL, y = "Methylation (%)") +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "02_Global_methylation_violin.png"),
       p_violin, width = 10, height = 6, dpi = 300)

# Density plot
p_density <- ggplot(long_meth, aes(x = methylation, colour = sample)) +
  geom_density(linewidth = 0.9) +
  scale_colour_manual(values = sample_colours) +
  theme_classic(base_size = 13) +
  labs(title = "Methylation Density Distribution",
       x = "Methylation (%)", y = "Density") +
  geom_vline(xintercept = c(20, 80), linetype = "dashed",
             colour = "grey50", linewidth = 0.5)

ggsave(file.path(out_dir, "03_Methylation_density.png"),
       p_density, width = 9, height = 5, dpi = 300)

# ─────────────────────────────────────────────────────────
#  STEP 5 — SAMPLE CORRELATION HEATMAP
# ─────────────────────────────────────────────────────────

cat("[5/16] Sample correlation heatmap...\n")

meth_mat <- as.matrix(merged[, .(C1, C2, C3, T1, T2, T3)])
cor_mat   <- cor(meth_mat, method = "pearson", use = "complete.obs")

ann_col <- data.frame(
  Group = sample_meta$group,
  row.names = sample_meta$sample_id
)

pheatmap(
  cor_mat,
  annotation_col    = ann_col,
  annotation_colors = list(Group = group_colours),
  color             = colorRampPalette(c("#053061","white","#67001F"))(100),
  breaks            = seq(0.8, 1, length.out = 101),
  display_numbers   = TRUE,
  number_format     = "%.3f",
  fontsize_number   = 9,
  main              = "Pearson Correlation — CpG Methylation",
  filename          = file.path(out_dir, "04_Sample_correlation_heatmap.png"),
  width = 7, height = 6
)

# ─────────────────────────────────────────────────────────
#  STEP 6 — PCA
# ─────────────────────────────────────────────────────────

cat("[6/16] Principal component analysis...\n")

# Remove zero-variance CpGs
meth_t <- t(meth_mat)
meth_t <- meth_t[, apply(meth_t, 2, var) > 0]

pca     <- prcomp(meth_t, scale. = TRUE)
var_exp <- round(100 * summary(pca)$importance[2, 1:4], 1)

pca_df <- data.frame(
  pca$x[, 1:4],
  sample  = sample_meta$sample_id,
  group   = sample_meta$group,
  patient = sample_meta$patient
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = group, shape = patient)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = sample), size = 3.5, show.legend = FALSE) +
  scale_colour_manual(values = group_colours) +
  theme_classic(base_size = 13) +
  labs(
    title    = "PCA — RRBS Methylation Profiles",
    subtitle = "GSE216050 | Prostate Cancer",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  )

ggsave(file.path(out_dir, "05_PCA_plot.png"),
       p_pca, width = 8, height = 6, dpi = 300)

# ─────────────────────────────────────────────────────────
#  STEP 7 — DIFFERENTIAL METHYLATION (t-test + FDR)
# ─────────────────────────────────────────────────────────

cat("[7/16] Differential methylation testing (t-test + BH FDR)...\n")

ctrl_cols <- c("C1","C2","C3")
trt_cols  <- c("T1","T2","T3")

merged[, control_mean := rowMeans(.SD), .SDcols = ctrl_cols]
merged[, treated_mean := rowMeans(.SD), .SDcols = trt_cols]
merged[, diff         := treated_mean - control_mean]

# Paired t-test (matched patients: PM154, PM12, PM21)
pvals <- apply(
  as.matrix(merged[, .(C1, C2, C3, T1, T2, T3)]),
  1,
  function(x) {
    ctrl <- x[1:3]
    trt  <- x[4:6]
    if (sd(c(ctrl, trt)) == 0) return(1)
    tryCatch(
      t.test(ctrl, trt, paired = TRUE)$p.value,
      error = function(e) 1
    )
  }
)

merged[, pvalue := pvals]
merged[, fdr    := p.adjust(pvalue, method = "BH")]
merged[, logP   := -log10(pvalue + 1e-10)]

# ─────────────────────────────────────────────────────────
#  STEP 8 — VOLCANO PLOT
# ─────────────────────────────────────────────────────────

cat("[8/16] Volcano plot...\n")

volcano_df <- as.data.frame(merged[, .(chr, start, diff, pvalue, fdr, logP)])
volcano_df$status <- "NS"
volcano_df$status[volcano_df$fdr    < FDR_CUTOFF &
                    abs(volcano_df$diff) > DIFF_THRESHOLD] <- "Significant"
volcano_df$status[volcano_df$diff   < -DIFF_THRESHOLD &
                    volcano_df$fdr    < FDR_CUTOFF]        <- "Hypomethylated"
volcano_df$status[volcano_df$diff   > DIFF_THRESHOLD  &
                    volcano_df$fdr    < FDR_CUTOFF]        <- "Hypermethylated"

status_cols <- c(
  NS              = "grey70",
  Significant     = "#FDAE61",
  Hypomethylated  = "#2166AC",
  Hypermethylated = "#D6604D"
)

n_hyper <- sum(volcano_df$status == "Hypermethylated")
n_hypo  <- sum(volcano_df$status == "Hypomethylated")

p_volcano <- ggplot(volcano_df,
                    aes(x = diff, y = logP, colour = status)) +
  geom_point(alpha = 0.35, size = 0.7) +
  scale_colour_manual(values = status_cols) +
  geom_vline(xintercept = c(-DIFF_THRESHOLD, DIFF_THRESHOLD),
             linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(PVAL_CUTOFF),
             linetype = "dashed", colour = "black", linewidth = 0.5) +
  annotate("text", x = -45, y = max(volcano_df$logP, na.rm=TRUE) * 0.95,
           label = paste0("Hypomethylated\nn = ", n_hypo),
           colour = "#2166AC", fontface = "bold", size = 4) +
  annotate("text", x =  45, y = max(volcano_df$logP, na.rm=TRUE) * 0.95,
           label = paste0("Hypermethylated\nn = ", n_hyper),
           colour = "#D6604D", fontface = "bold", size = 4) +
  theme_classic(base_size = 13) +
  labs(
    title    = "Differential CpG Methylation — Decitabine vs Control",
    subtitle = paste0("FDR < ", FDR_CUTOFF, "  |  |Δ Methylation| > ", DIFF_THRESHOLD, "%"),
    x        = "Methylation Difference (Treated − Control, %)",
    y        = expression(-log[10](p~value)),
    colour   = NULL
  ) +
  theme(legend.position = "top")

ggsave(file.path(out_dir, "06_Volcano_plot.png"),
       p_volcano, width = 10, height = 7, dpi = 300)

# ─────────────────────────────────────────────────────────
#  STEP 9 — SIGNIFICANT CpGs
# ─────────────────────────────────────────────────────────

cat("[9/16] Extracting significant CpGs...\n")

sig_cpg <- merged[fdr < FDR_CUTOFF & abs(diff) > DIFF_THRESHOLD]
cat("   Significant CpGs (FDR<", FDR_CUTOFF, ", |Δ|>", DIFF_THRESHOLD, "%):",
    nrow(sig_cpg), "\n")
cat("   Hypermethylated (treated):", sum(sig_cpg$diff > 0), "\n")
cat("   Hypomethylated  (treated):", sum(sig_cpg$diff < 0), "\n")

write.csv(
  sig_cpg[order(fdr)],
  file.path(out_dir, "07_Significant_CpGs.csv"),
  row.names = FALSE
)

# ─────────────────────────────────────────────────────────
#  STEP 10 — DMR CALLING (sliding-window approach)
# ─────────────────────────────────────────────────────────

cat("[10/16] Calling DMRs by sliding window...\n")

sig_gr <- GRanges(
  seqnames = sig_cpg$chr,
  ranges   = IRanges(start = sig_cpg$start, end = sig_cpg$end),
  diff     = sig_cpg$diff,
  fdr      = sig_cpg$fdr
)

# Reduce nearby CpGs into windows
sig_gr_reduced <- reduce(resize(sig_gr, width = DMR_WIN_SIZE, fix = "center"))

# Annotate each window with CpG count + mean diff
hits <- findOverlaps(sig_gr_reduced, sig_gr)

# Build stats per window safely
query_idx   <- queryHits(hits)
subject_idx <- subjectHits(hits)

window_n_cpgs    <- as.integer(tapply(subject_idx, query_idx, length))
window_mean_diff <- as.numeric(tapply(sig_gr$diff[subject_idx], query_idx, mean))
window_min_fdr   <- as.numeric(tapply(sig_gr$fdr[subject_idx],  query_idx, min))

dmr_dt <- as.data.table(sig_gr_reduced)
dmr_dt[, n_cpgs    := window_n_cpgs]
dmr_dt[, mean_diff := round(window_mean_diff, 3)]
dmr_dt[, min_fdr   := signif(window_min_fdr, 4)]

# Filter to windows with enough CpGs
dmr_dt <- dmr_dt[n_cpgs >= DMR_MIN_CPGS]
setnames(dmr_dt, "seqnames", "chr")

cat("   DMRs identified (>=", DMR_MIN_CPGS, "CpGs per window):", nrow(dmr_dt), "\n")

write.csv(dmr_dt[order(min_fdr)],
          file.path(out_dir, "08_DMR_list.csv"),
          row.names = FALSE)

# ─────────────────────────────────────────────────────────
#  STEP 11 — GENOMIC ANNOTATION
# ─────────────────────────────────────────────────────────

cat("[11/16] Annotating significant CpGs to genomic features...\n")

txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)

# Promoters (2kb upstream / 500bp downstream of TSS)
proms   <- promoters(txdb, upstream = 2000, downstream = 500)
exons   <- exons(txdb)
introns <- gaps(reduce(exons))

sig_gr2 <- GRanges(
  seqnames = sig_cpg$chr,
  ranges   = IRanges(start = sig_cpg$start, end = sig_cpg$end),
  diff     = sig_cpg$diff
)

annotate_feature <- function(query, feature, label) {
  ov <- overlapsAny(query, feature)
  sum(ov)
}

annot_counts <- data.frame(
  Feature = c("Promoter (±2kb TSS)", "Exon", "Intron", "Intergenic"),
  Count   = c(
    annotate_feature(sig_gr2, proms,   "Promoter"),
    annotate_feature(sig_gr2, exons,   "Exon"),
    annotate_feature(sig_gr2, introns, "Intron"),
    nrow(sig_cpg) - annotate_feature(sig_gr2, c(proms, exons, introns), "all")
  )
)
annot_counts$Percent <- round(100 * annot_counts$Count / nrow(sig_cpg), 1)

p_annot <- ggplot(annot_counts,
                  aes(x = reorder(Feature, -Count),
                      y = Count, fill = Feature)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(Count, "\n(", Percent, "%)")),
            vjust = -0.3, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 13) +
  labs(
    title = "Genomic Distribution of Significant CpGs",
    x     = NULL,
    y     = "Number of CpGs"
  )

ggsave(file.path(out_dir, "09_Genomic_annotation_barplot.png"),
       p_annot, width = 7, height = 5, dpi = 300)

# Promoter-specific DMR genes
hits_gene  <- findOverlaps(sig_gr2, genes)
gene_ids   <- unique(genes$gene_id[subjectHits(hits_gene)])
gene_syms  <- bitr(gene_ids, fromType = "ENTREZID",
                   toType = c("SYMBOL","GENENAME"),
                   OrgDb = org.Hs.eg.db)

# Annotate direction
hits2         <- findOverlaps(sig_gr2, genes)
gene_id_full  <- genes$gene_id[subjectHits(hits2)]
diff_vals     <- sig_gr2$diff[queryHits(hits2)]
gene_dir      <- data.table(entrez = gene_id_full, diff = diff_vals)
gene_dir_agg  <- gene_dir[, .(mean_diff = mean(diff)), by = entrez]
gene_syms     <- merge(gene_syms,
                       gene_dir_agg,
                       by.x = "ENTREZID",
                       by.y = "entrez",
                       all.x = TRUE)

write.csv(gene_syms,
          file.path(out_dir, "10_DMR_gene_list.csv"),
          row.names = FALSE)

cat("   Genes overlapping DMRs:", nrow(gene_syms), "\n")

# ─────────────────────────────────────────────────────────
#  STEP 12 — TOP CpG HEATMAP
# ─────────────────────────────────────────────────────────

cat("[12/16] Heatmap of top 100 differentially methylated CpGs...\n")

top100 <- merged[order(fdr)][1:min(100, nrow(merged))]

heat_mat <- as.matrix(top100[, .(C1, C2, C3, T1, T2, T3)])
rownames(heat_mat) <- paste0(top100$chr, ":", top100$start)

ann_col_heat <- data.frame(
  Group   = sample_meta$group,
  Patient = sample_meta$patient,
  row.names = sample_meta$sample_id
)
ann_colors <- list(
  Group   = group_colours,
  Patient = c(PM154 = "#7B2D8B", PM12 = "#1A9850", PM21 = "#D73027")
)

pheatmap(
  heat_mat,
  annotation_col    = ann_col_heat,
  annotation_colors = ann_colors,
  scale             = "row",
  clustering_method = "ward.D2",
  color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  show_rownames     = FALSE,
  fontsize_col      = 10,
  main              = "Top 100 Differentially Methylated CpGs\n(scaled by row)",
  filename          = file.path(out_dir, "11_Top100_CpG_heatmap.png"),
  width = 8, height = 8
)

# ─────────────────────────────────────────────────────────
#  STEP 13 — B7-H3 (CD276) LOCUS SPOTLIGHT
#   This is the key gene in the study — decitabine upregulates
#   B7-H3 via CpG island demethylation
# ─────────────────────────────────────────────────────────

cat("[13/16] B7-H3 (CD276) locus-specific methylation...\n")

# CD276 / B7-H3 is on chr15: 38,534,000 – 38,630,000 (hg38)
b7h3_range <- GRanges("chr15", IRanges(38534000, 38630000))

all_gr <- GRanges(
  seqnames = merged$chr,
  ranges   = IRanges(start = merged$start, end = merged$end),
  C1 = merged$C1, C2 = merged$C2, C3 = merged$C3,
  T1 = merged$T1, T2 = merged$T2, T3 = merged$T3,
  diff = merged$diff
)

b7h3_hits <- subsetByOverlaps(all_gr, b7h3_range)

if (length(b7h3_hits) > 0) {
  b7h3_df <- as.data.table(b7h3_hits)[, .(start, C1, C2, C3, T1, T2, T3, diff)]
  b7h3_long <- melt(b7h3_df, id.vars = c("start","diff"),
                    variable.name = "sample", value.name = "methylation")
  b7h3_long[, group := ifelse(sample %in% c("C1","C2","C3"),
                              "Control","Decitabine")]
  
  p_b7h3 <- ggplot(b7h3_long,
                   aes(x = start, y = methylation,
                       colour = group, group = sample)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = group_colours) +
    theme_classic(base_size = 12) +
    labs(
      title    = "B7-H3 (CD276) Locus — CpG Methylation Profile",
      subtitle = "chr15:38,534,000–38,630,000 | Decitabine effect on B7-H3 promoter",
      x        = "Genomic Position (chr15)",
      y        = "Methylation (%)",
      colour   = NULL
    ) +
    scale_x_continuous(labels = comma)
  
  ggsave(file.path(out_dir, "12_B7H3_locus_methylation.png"),
         p_b7h3, width = 10, height = 5, dpi = 300)
} else {
  cat("   No CpGs found at the B7-H3 locus in common set — check coordinates.\n")
}

# ─────────────────────────────────────────────────────────
#  STEP 14 — GO ENRICHMENT (Biological Process)
# ─────────────────────────────────────────────────────────

cat("[14/16] GO enrichment analysis...\n")

gene_list_entrez <- gene_syms$ENTREZID

ego <- enrichGO(
  gene          = gene_list_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.10,
  readable      = TRUE
)

if (!is.null(ego) && nrow(ego) > 0) {
  p_go <- dotplot(ego, showCategory = 20, font.size = 10) +
    labs(title = "GO Biological Process Enrichment",
         subtitle = "Genes overlapping significant DMRs") +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(out_dir, "13_GO_enrichment_dotplot.png"),
         p_go, width = 10, height = 9, dpi = 300)
  
  write.csv(as.data.frame(ego),
            file.path(out_dir, "14_GO_enrichment_table.csv"),
            row.names = FALSE)
  
  cat("   GO terms enriched:", nrow(ego), "\n")
} else {
  cat("   No significant GO terms found — try loosening FDR threshold.\n")
}

# ─────────────────────────────────────────────────────────
#  STEP 15 — KEGG PATHWAY ENRICHMENT
# ─────────────────────────────────────────────────────────

cat("[15/16] KEGG pathway enrichment...\n")

kegg <- enrichKEGG(
  gene          = gene_list_entrez,
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.10
)

if (!is.null(kegg) && nrow(kegg) > 0) {
  p_kegg <- dotplot(kegg, showCategory = 20, font.size = 10) +
    labs(title = "KEGG Pathway Enrichment",
         subtitle = "Genes overlapping significant DMRs") +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(out_dir, "15_KEGG_enrichment_dotplot.png"),
         p_kegg, width = 10, height = 9, dpi = 300)
  
  write.csv(as.data.frame(kegg),
            file.path(out_dir, "16_KEGG_enrichment_table.csv"),
            row.names = FALSE)
  
  cat("   KEGG pathways enriched:", nrow(kegg), "\n")
} else {
  cat("   No significant KEGG pathways found.\n")
}

# ─────────────────────────────────────────────────────────
#  STEP 16 — FINAL SUMMARY REPORT
# ─────────────────────────────────────────────────────────

cat("[16/16] Writing summary report...\n")

summary_report <- data.frame(
  Metric = c(
    "Dataset",
    "Common CpGs (all 6 samples)",
    "Coverage threshold",
    "Differential test",
    "FDR method",
    "Sig. CpGs (FDR<0.05, |Δ|>20%)",
    "Hypermethylated in Decitabine",
    "Hypomethylated in Decitabine",
    "DMRs identified",
    "Genes overlapping DMRs",
    "GO BP terms enriched",
    "KEGG pathways enriched"
  ),
  Value = c(
    "GSE216050",
    format(nrow(merged), big.mark=","),
    paste0("≥", MIN_COVERAGE, " reads"),
    "Paired t-test",
    "Benjamini-Hochberg",
    format(nrow(sig_cpg), big.mark=","),
    format(sum(sig_cpg$diff > 0), big.mark=","),
    format(sum(sig_cpg$diff < 0), big.mark=","),
    format(nrow(dmr_dt), big.mark=","),
    format(nrow(gene_syms), big.mark=","),
    ifelse(!is.null(ego)  && nrow(ego)  > 0, nrow(ego),  "0"),
    ifelse(!is.null(kegg) && nrow(kegg) > 0, nrow(kegg), "0")
  )
)

write.csv(summary_report,
          file.path(out_dir, "00_Analysis_summary.csv"),
          row.names = FALSE)

cat("\n============================================================\n")
cat(" ANALYSIS COMPLETE\n")
cat(" All outputs saved to:", out_dir, "\n")
cat("============================================================\n")
print(summary_report)

# List all output files
cat("\nFiles generated:\n")
out_files <- list.files(out_dir, full.names = FALSE)
for (f in out_files) cat("  →", f, "\n")
