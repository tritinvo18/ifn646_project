# Task 1: RNA-seq Differential Expression Analysis & Tool Consensus
# Dataset: GSE159717 (SARS-CoV-2 infection of human islets, 5 dpi)
# Tools compared: DESeq2 vs edgeR vs limma-voom

# Suppress package startup messages
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(grid)
  library(ashr)
  library(limma)
})

# Set working directory to project root if needed
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
# Create output directories if they don't exist
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Load Data & Experimental Metadata
# ------------------------------------------------------------------------------
counts_file <- "gse159717_rnaseq_deseq_5dpi_counts_raw.tsv"
cat("Loading raw counts from:", counts_file, "\n")
raw_data <- read.delim(counts_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")

cat("Raw count matrix dimensions:", nrow(raw_data), "genes by", ncol(raw_data), "columns\n")

# Extract gene annotation columns
gene_annot <- raw_data[, c("gene_id", "gene_version", "gene_name", "gene_biotype", "gene_source", "location")]
rownames(gene_annot) <- gene_annot$gene_id

# Extract chromosome number from location string (e.g. chr1:11869-14409:+ -> 1)
gene_annot$chromosome <- gsub("^chr([^:]+):.*$", "\\1", gene_annot$location)

# Define sample columns
sample_cols <- c("S_2_mock_5dpi_S70002", "S_2_SARS_5dpi_S70003", "S_2_Rem_5dpi_S70001",
                 "S_3_mock_5dpi_S69997", "S_3_SARS_5dpi_S69996", "S_3_Rem_5dpi_S69995")

# Construct metadata dataframe
col_data <- data.frame(
  sample_id = sample_cols,
  donor     = factor(c("Donor_2", "Donor_2", "Donor_2", "Donor_3", "Donor_3", "Donor_3")),
  condition = factor(c("mock", "SARS", "Rem", "mock", "SARS", "Rem"), levels = c("mock", "SARS", "Rem")),
  row.names = sample_cols,
  stringsAsFactors = FALSE
)
cat("\nExperimental Design Metadata:\n")
print(col_data)

# Raw count matrix
count_mat <- as.matrix(raw_data[, sample_cols])
rownames(count_mat) <- raw_data$gene_id


# EDA - Data Quality Assessment
# ------------------------------------------------------------------------------

# Filter genes: retain genes with >= 10 counts in at least 2 samples
keep_genes <- rowSums(count_mat >= 10) >= 2
filtered_counts <- count_mat[keep_genes, ]
filtered_annot  <- gene_annot[keep_genes, ]

cat("\nFiltering summary:\n")
cat("Total genes before filtering:", nrow(count_mat), "\n")
cat("Genes retained after filtering (>= 10 counts in >= 2 samples):", nrow(filtered_counts), "\n")
cat("Genes filtered out:", nrow(count_mat) - nrow(filtered_counts), "\n\n")

# DESeq2 Pipeline
# ------------------------------------------------------------------------------

cat("--- DESeq2 ---\n")
# Paired design: ~ donor + condition (controls for biological variation across donors)
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData   = col_data,
  design    = ~ donor + condition
)

# Run DESeq pipeline (Size factor estimation, dispersion estimation, Wald test)
dds <- DESeq(dds, quiet = FALSE)

# Variance-stabilizing transformation for EDA
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

# PCA Plot
pca_data <- plotPCA(rld, intgroup = c("condition", "donor"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = donor)) +
  geom_point(size = 5, stroke = 1.2) +
  scale_color_manual(values = c("mock" = "#2b83ba", "SARS" = "#d7191c", "Rem" = "#fdae61")) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("Principal Component Analysis (PCA) - GSE159717") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/01_pca_plot.png", pca_plot, width = 7, height = 5, dpi = 300)
cat("Saved PCA plot to results/figures/01_pca_plot.png\n")

# Sample Correlation Heatmap
sample_dists <- dist(t(assay(rld)))
sample_dist_mat <- as.matrix(sample_dists)
rownames(sample_dist_mat) <- paste(col_data$donor, col_data$condition, sep = "_")
colnames(sample_dist_mat) <- paste(col_data$donor, col_data$condition, sep = "_")
heat_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png("results/figures/02_sample_distance_heatmap.png", width = 1800, height = 1500, res = 300)
pheatmap(sample_dist_mat,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = heat_colors,
         main = "Sample-to-Sample Distance Heatmap (rlog)")
dev.off()
cat("Saved sample distance heatmap to results/figures/02_sample_distance_heatmap.png\n")

# Sample correlation heatmap
cor_mat <- cor(assay(rld), method = "pearson")
rownames(cor_mat) <- paste(col_data$donor, col_data$condition, sep = "_")
colnames(cor_mat) <- paste(col_data$donor, col_data$condition, sep = "_")
png("results/figures/02b_sample_correlation_heatmap.png", width = 1800, height = 1500, res = 300)
pheatmap(cor_mat,
         col = colorRampPalette(brewer.pal(9, "YlOrRd"))(255),
         main = "Pairwise Sample Correlation Matrix (rlog)")
dev.off()
cat("Saved correlation heatmap to results/figures/02b_sample_correlation_heatmap.png\n")

# Dispersion Estimates Plot
png("results/figures/03_deseq2_dispersion_estimates.png", width = 2100, height = 1800, res = 300)
plotDispEsts(dds, main = "DESeq2 Dispersion Estimates Shrinkage")
dev.off()
cat("Saved dispersion plot to results/figures/03_deseq2_dispersion_estimates.png\n")

# Extract DESeq2 Results (SARS vs mock)
res_deseq_raw <- results(dds, contrast = c("condition", "SARS", "mock"), alpha = 0.05)

# Perform Log2 Fold Change Shrinkage (ashr)
res_deseq_shrunk <- lfcShrink(dds, contrast = c("condition", "SARS", "mock"),
                              res = res_deseq_raw, type = "ashr")

# MA Plots: Unshrunken vs Shrunken
png("results/figures/04_deseq2_ma_plots.png", width = 3000, height = 1500, res = 300)
par(mfrow = c(1, 2))
plotMA(res_deseq_raw, ylim = c(-6, 6), main = "DESeq2 MA Plot (Unshrunken LFC)")
plotMA(res_deseq_shrunk, ylim = c(-6, 6), main = "DESeq2 MA Plot (ashr Shrunken LFC)")
dev.off()
cat("Saved MA plots to results/figures/04_deseq2_ma_plots.png\n")

# Convert to data frame
deseq_df <- as.data.frame(res_deseq_shrunk)
deseq_df$gene_id <- rownames(deseq_df)
deseq_df <- left_join(deseq_df, filtered_annot, by = "gene_id")

# Save DESeq2 full results
write.csv(deseq_df, "results/tables/deseq2_sars_vs_mock_all_genes.csv", row.names = FALSE)
cat("Saved DESeq2 results table to results/tables/deseq2_sars_vs_mock_all_genes.csv\n")


# edgeR Pipeline
# ------------------------------------------------------------------------------
cat("\n--- edgeR (GLM Likelihood Ratio Test) ---\n")
# Create DGEList
dge <- DGEList(counts = filtered_counts, genes = filtered_annot)

# TMM Normalization
dge <- normLibSizes(dge) # or calcNormFactors
cat("TMM Normalization Factors:\n")
print(dge$samples)

# Model design matrix
design_mat <- model.matrix(~ donor + condition, data = col_data)
cat("edgeR Design Matrix:\n")
print(design_mat)

# Estimate Dispersions (common, trended, tagwise)
dge <- estimateDisp(dge, design = design_mat)

# Plot Biological Coefficient of Variation (BCV)
png("results/figures/05_edger_bcv_plot.png", width = 2100, height = 1800, res = 300)
plotBCV(dge, main = "edgeR Biological Coefficient of Variation (BCV)")
dev.off()
cat("Saved edgeR BCV plot to results/figures/05_edger_bcv_plot.png\n")

# Fit negative binomial GLM
fit_edger <- glmFit(dge, design = design_mat)

# Likelihood Ratio Test for condition SARS (coef = "conditionSARS")
lrt_sars <- glmLRT(fit_edger, coef = "conditionSARS")
edger_top <- topTags(lrt_sars, n = Inf)$table

# edgeR MA/MD plot
png("results/figures/06_edger_md_plot.png", width = 2100, height = 1800, res = 300)
plotMD(lrt_sars, status = decideTests(lrt_sars, p.value = 0.05, lfc = 1),
       main = "edgeR Mean-Difference (MD) Plot (SARS vs mock)")
dev.off()
cat("Saved edgeR MD plot to results/figures/06_edger_md_plot.png\n")

# Save edgeR full results
write.csv(edger_top, "results/tables/edger_sars_vs_mock_all_genes.csv", row.names = FALSE)
cat("Saved edgeR results table to results/tables/edger_sars_vs_mock_all_genes.csv\n")


# limma-voom Pipeline
# ------------------------------------------------------------------------------

cat("\n--- limma-voom ---\n")
v <- voom(dge, design_mat, plot=FALSE)

png("results/figures/07_limma_voom_plot.png", width = 2100, height = 1800, res = 300)
voom(dge, design_mat, plot=TRUE)
dev.off()
cat("Saved limma-voom trend plot to results/figures/07_limma_voom_plot.png\n")

fit_limma <- lmFit(v, design_mat)
fit_limma <- eBayes(fit_limma)

limma_res <- topTable(fit_limma, coef="conditionSARS", n=Inf)
limma_res$gene_id <- rownames(limma_res)

write.csv(limma_res, "results/tables/limma_voom_sars_vs_mock_all_genes.csv", row.names = FALSE)
cat("Saved limma-voom results table to results/tables/limma_voom_sars_vs_mock_all_genes.csv\n")


# Consensus Analysis & Tool Comparison
# ------------------------------------------------------------------------------
cat("\n--- Conducting Consensus Analysis (DESeq2 vs edgeR) ---\n")

# Criteria thresholds:
# 1. Adjusted p-value (FDR) < 0.05
# 2. Absolute Log2 Fold Change >= 1.0 (2-fold change)
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0

# Identify significant genes in DESeq2
deseq_sig <- deseq_df %>%
  filter(!is.na(padj) & padj < PADJ_CUTOFF & abs(log2FoldChange) >= LFC_CUTOFF)

# Identify significant genes in edgeR
edger_sig <- edger_top %>%
  filter(!is.na(FDR) & FDR < PADJ_CUTOFF & abs(logFC) >= LFC_CUTOFF)

# Identify significant genes in limma
limma_sig <- limma_res %>%
  filter(!is.na(adj.P.Val) & adj.P.Val < PADJ_CUTOFF & abs(logFC) >= LFC_CUTOFF)

cat("DESeq2 significant genes:", nrow(deseq_sig), "\n")
cat("edgeR  significant genes:", nrow(edger_sig), "\n")
cat("limma-voom significant genes:", nrow(limma_sig), "\n")

# Intersect significant gene IDs
common_ids_de <- intersect(deseq_sig$gene_id, edger_sig$gene_id)
common_ids_dl <- intersect(deseq_sig$gene_id, limma_sig$gene_id)
common_ids_el <- intersect(edger_sig$gene_id, limma_sig$gene_id)

# limma-voom returned 0 genes, so a 3-way intersection would be empty.
# We will define the strict consensus based on the best pair (DESeq2 and edgeR).
common_ids <- common_ids_de
cat("Strict Consensus Genes (DESeq2 & edgeR):", length(common_ids), "\n")

# Evaluate slightly more relaxed threshold (FDR < 0.10, |log2FC| >= 0.58)
deseq_relaxed <- deseq_df %>%
  filter(!is.na(padj) & padj < 0.10 & abs(log2FoldChange) >= 0.58)
edger_relaxed <- edger_top %>%
  filter(!is.na(FDR) & FDR < 0.10 & abs(logFC) >= 0.58)
relaxed_common_ids <- intersect(deseq_relaxed$gene_id, edger_relaxed$gene_id)
cat("Relaxed Consensus Genes (FDR < 0.10, |log2FC| >= 0.58):", length(relaxed_common_ids), "\n")

# Log2 Fold Change Correlation Plot
# Merge all tested genes
comparison_df <- inner_join(
  deseq_df %>% dplyr::select(gene_id, gene_name, deseq_lfc = log2FoldChange, deseq_padj = padj, chromosome, gene_biotype),
  edger_top %>% dplyr::select(gene_id, edger_lfc = logFC, edger_fdr = FDR),
  by = "gene_id"
) %>% inner_join(
  limma_res %>% dplyr::select(gene_id, limma_lfc = logFC, limma_padj = adj.P.Val),
  by = "gene_id"
)

comparison_df$status <- "Not DE"
comparison_df$status[comparison_df$gene_id %in% deseq_sig$gene_id] <- "DESeq2 only"
comparison_df$status[comparison_df$gene_id %in% edger_sig$gene_id] <- "edgeR only"
comparison_df$status[comparison_df$gene_id %in% common_ids]       <- "Consensus DE"
comparison_df$status <- factor(comparison_df$status, levels = c("Not DE", "DESeq2 only", "edgeR only", "Consensus DE"))

# Compute Spearman rank correlation of log2FC
cor_lfc <- cor(comparison_df$deseq_lfc, comparison_df$edger_lfc, method = "spearman", use = "complete.obs")
cat("Spearman correlation of Log2 Fold Changes between DESeq2 and edgeR:", round(cor_lfc, 4), "\n")

lfc_scatter <- ggplot(comparison_df, aes(x = deseq_lfc, y = edger_lfc, color = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = c("Not DE" = "#bdbdbd", "DESeq2 only" = "#3182bd", "edgeR only" = "#de2d26", "Consensus DE" = "#31a354")) +
  geom_text_repel(data = comparison_df %>% filter(status == "Consensus DE" | gene_name %in% c("ISG15", "IFI6", "GBP4", "PLEKHM3")),
                  aes(label = gene_name), size = 3.5, max.overlaps = 20, color = "black", box.padding = 0.5) +
  labs(
    x = "DESeq2 Log2 Fold Change (shrunken)",
    y = "edgeR Log2 Fold Change",
    title = "Cross-Tool Concordance of Estimated Effect Sizes",
    subtitle = paste0("Spearman Rank Correlation r = ", round(cor_lfc, 3))
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("results/figures/08_lfc_concordance_scatter.png", lfc_scatter, width = 8, height = 7, dpi = 300)
cat("Saved LFC concordance scatter plot to results/figures/08_lfc_concordance_scatter.png\n")

# Volcano Plots with Highlights
volcano_df <- deseq_df %>%
  mutate(
    neg_log10_p = -log10(padj),
    significance = case_when(
      gene_id %in% common_ids ~ "Consensus DE",
      padj < PADJ_CUTOFF & abs(log2FoldChange) >= LFC_CUTOFF ~ "DESeq2 DE",
      TRUE ~ "NS"
    )
  )

volcano_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_p, color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("NS" = "#bdbdbd", "DESeq2 DE" = "#3182bd", "Consensus DE" = "#31a354")) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = volcano_df %>% filter(significance == "Consensus DE" | gene_name %in% c("ISG15", "IFI6", "GBP4", "PLEKHM3")),
                  aes(label = gene_name), size = 3.8, max.overlaps = 25, color = "black", box.padding = 0.5) +
  labs(
    x = "Log2 Fold Change (SARS vs mock)",
    y = "-Log10 Adjusted P-Value",
    title = "DESeq2 Volcano Plot (SARS vs mock)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("results/figures/09_volcano_plot.png", volcano_plot, width = 8, height = 7, dpi = 300)
cat("Saved Volcano plot to results/figures/09_volcano_plot.png\n")


# Build Comprehensive Consensus Tables & Chromosome 1 / 2 Subsets
# ------------------------------------------------------------------------------
cat("\n--- Consensus Tables ---\n")

# Merge metrics for consensus genes
consensus_master <- comparison_df %>%
  filter(status == "Consensus DE") %>%
  left_join(filtered_annot %>% dplyr::select(gene_id, location), by = "gene_id") %>%
  arrange(deseq_padj)

# Save strict consensus table
write.csv(consensus_master, "results/tables/consensus_de_genes_strict.csv", row.names = FALSE)
cat("Saved strict consensus genes to results/tables/consensus_de_genes_strict.csv\n")
print(consensus_master)

# List Chromosome 1 and Chromosome 2 prioritized table for Tasks 2 and 3!
# includes both strict consensus and top high-significance candidates on Chr 1 and Chr 2
chr1_2_table <- comparison_df %>%
  filter(chromosome %in% c("1", "2") & gene_biotype == "protein_coding") %>%
  filter(deseq_padj < 0.10 | edger_fdr < 0.10) %>%
  left_join(filtered_annot %>% dplyr::select(gene_id, location), by = "gene_id") %>%
  arrange(deseq_padj)

# Add suggested CRISPR modality
chr1_2_table$crispr_modality <- ifelse(chr1_2_table$deseq_lfc > 0, "CRISPRi (Repression)", "CRISPRa (Activation)")

write.csv(chr1_2_table, "results/tables/target_candidates_chr1_chr2.csv", row.names = FALSE)
cat("\nSaved prioritized Chr1 / Chr2 target candidate table to results/tables/target_candidates_chr1_chr2.csv\n")
cat("\nTop Candidate Genes on Chromosome 1 and Chromosome 2:\n")
print(head(chr1_2_table, 15))


# Heatmap of Top Differentially Expressed Genes Across All Conditions
# ------------------------------------------------------------------------------
cat("\n Expression heatmap for top candidates \n")
top_genes_for_heatmap <- unique(c(consensus_master$gene_id, head(chr1_2_table$gene_id, 15)))

# Extract normalized counts
norm_counts <- counts(dds, normalized = TRUE)[top_genes_for_heatmap, ]
# Z-score scale by row
norm_counts_scaled <- t(scale(t(norm_counts)))

# Match row labels with gene symbols
rownames(norm_counts_scaled) <- gene_annot[top_genes_for_heatmap, "gene_name"]
colnames(norm_counts_scaled) <- paste(col_data$donor, col_data$condition, sep = "_")

# Annotation for columns
annotation_col <- data.frame(
  Condition = col_data$condition,
  Donor     = col_data$donor,
  row.names = colnames(norm_counts_scaled)
)
annotation_colors <- list(
  Condition = c("mock" = "#2b83ba", "SARS" = "#d7191c", "Rem" = "#fdae61"),
  Donor     = c("Donor_2" = "#7570b3", "Donor_3" = "#e7298a")
)

png("results/figures/10_top_genes_expression_heatmap.png", width = 2400, height = 2400, res = 300)
pheatmap(norm_counts_scaled,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = TRUE,
         main = "Expression Patterns of Top Candidate Genes (Z-score Normalized)")
dev.off()
cat("Saved expression heatmap to results/figures/10_top_genes_expression_heatmap.png\n")
