## Setup
### Bioconductor and CRAN libraries used
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(edgeR)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")

####---- Limma Voom ----
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
data0 <- read.table("gse159717_rnaseq_deseq_5dpi_counts_raw.tsv", header=T, row.names=1) 
data <- data0[c("S_2_Rem_5dpi_S70001", "S_2_SARS_5dpi_S70003","S_2_mock_5dpi_S70002","S_3_Rem_5dpi_S69995","S_3_SARS_5dpi_S69996","S_3_mock_5dpi_S69997")]
d0 <- DGEList(data)
d0 <- normLibSizes(d0)
d0
dim(d0)

# filter low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
head(d)

snames <- colnames(data) # Sample names
grps <- substr(snames, 5, 8)
batch <- substr(snames, 3, 3)
grps

plotMDS(d)

mm <- model.matrix(~0 + grps) # + batch)
y <- voom(d, mm, plot = T)

# Fit a linear model
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(grpsSARS - grpsmock, levels = colnames(coef(fit)))
contr

# Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

# Empirical Bayes smoothing of standard errors (shrinks standard errors that are 
# much larger or smaller than those from other genes towards the average standard error) 
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

####---- EdgeR ----
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
data0 <- read.table("gse159717_rnaseq_deseq_5dpi_counts_raw.tsv", header=T, row.names=1) 
data <- data0[c("S_2_Rem_5dpi_S70001", "S_2_SARS_5dpi_S70003","S_2_mock_5dpi_S70002","S_3_Rem_5dpi_S69995","S_3_SARS_5dpi_S69996","S_3_mock_5dpi_S69997")]
# d0 <- DGEList(data)
# d0 <- normLibSizes(d0)
# d0
# dim(d0)
# 
# # filter low-expressed genes
# # cutoff <- 1
# # drop <- which(apply(cpm(d0), 1, max) < cutoff)
# d <- d0 #[-drop,] 
# dim(d) # number of genes left
# head(d)

# Groups and batchs
snames <- colnames(data) # Sample names
grps <- substr(snames, 5, 8)
batch <- substr(snames, 3, 3)
grps

d <- DGEList(data, group = grps)
d$samples

# Filter out lowly expressed genes
keep <- filterByExpr(d, group = grps)
y <- d[keep, , keep.lib.sizes=FALSE]
dim(y)

design <- model.matrix(~0 + grps)
design

y <- normLibSizes(y)
y$samples

y <- estimateDisp(y, design)
fit <- glmQLFit(y,design)

chosen_contrast <- makeContrasts(grpsSARS-grpsmock, levels=design)
qlf <- glmQLFTest(fit, contrast=chosen_contrast)
topTags(qlf, n = 20)


####---- DESeq ----

data0 <- read.table("gse159717_rnaseq_deseq_5dpi_counts_raw.tsv", header=T, row.names=1) 
data <- data0[c("S_2_Rem_5dpi_S70001", "S_2_SARS_5dpi_S70003","S_2_mock_5dpi_S70002","S_3_Rem_5dpi_S69995","S_3_SARS_5dpi_S69996","S_3_mock_5dpi_S69997")]

meta <- read.csv("RNA metadata.csv", header=T, row.names=1) 
meta$batch <- factor(meta$batch)

### Check that sample names match in both files
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ batch + condition)
dds <- estimateSizeFactors(dds)

### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

### Plot PCA 
plotPCA(rld, intgroup="condition")
plotPCA(rld, intgroup="batch") # concerning..
  
### Plot correlation
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
rld_cor <- cor(rld_mat)    ## cor() is a base R function
# Plot heatmap
pheatmap(rld_cor)
  
# Run base method from https://bioconductor.org/packages//release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
dds <- DESeq(dds)

# Plot dispersion
plotDispEsts(dds)

l_contrast = c("condition", "CovOnly", "Control")
res <- results(dds, contrast = l_contrast, alpha = 0.05)
res_shrunken <- lfcShrink(dds,contrast = l_contrast, res = res, type = "ashr")

resOrdered <- res_shrunken[order(res_shrunken$pvalue),]
summary(resOrdered)
resOrdered

# Plot log of fold changes
plotMA(resOrdered, ylim=c(-4,4))

plotCounts(dds,gene = "ENSG00000111335", intgroup = "condition" )

colnames(resOrdered)


###---- Compare Results ----

limma_voom_poss <- top.table %>%
  mutate(rn = row_number()) %>%
  filter(P.Value < 0.01) %>%
  select("rn", "P.Value", "adj.P.Val")
colnames(limma_voom_poss) <- c("lv_rn", "lv_pval", "lv_adj_pval")

edger_poss <- as.data.frame(topTags(qlf, n = 100)) %>%
  arrange(PValue) %>%
  mutate(rn = row_number()) %>%
  filter(PValue < 0.05) %>%
  select("rn", "PValue", "FDR")
colnames(edger_poss) <- c("ed_rn", "ed_pval", "ed_adj_pval")

deseq_poss <- as.data.frame(resOrdered) %>% #used Copilot for this syntax 
  mutate(rn = row_number()) %>%
  filter(padj < 0.05) %>%
  select("rn", "pvalue", "padj") 
colnames(deseq_poss) <- c("ds_rn", "ds_pval", "ds_adj_pval")


merged_df <- merge(
  deseq_poss, limma_voom_poss, by = "row.names", all = FALSE
  )
row.names(merged_df) <- merged_df$Row.names
merged_df <- merged_df[2:7]

merged_df <- merge(
    merged_df, edger_poss, by = "row.names", all = FALSE
  )  %>%
  mutate(avg_rn = (ds_rn+lv_rn+ed_rn)/3) %>%
  arrange(avg_rn)

merged_df
