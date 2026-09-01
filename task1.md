# Task 1: RNA-seq Differential Expression Analysis & Tool Consensus

## Overview
This report summarizes the differential expression analysis for Task 1 based on dataset **GSE159717** (SARS-CoV-2 infection of human islets at 5 days post-infection). The goal was to identify differentially expressed genes (DEGs) between SARS-CoV-2 infected samples and mock controls by running two independent analytical pipelines and evaluating their consensus.

## Methods
Three major RNA-seq analysis pipelines were employed, accounting for both `condition` (SARS vs mock) and `donor` (to control for biological variation). DESeq2, edgeR, and limma-voom are all highly suitable for this problem. DESeq2 and edgeR model raw count data using a negative binomial distribution to capture biological overdispersion, while limma-voom transforms counts to log2-counts per million (logCPM) and estimates precision weights based on the mean-variance trend. All tools implement sophisticated empirical Bayes techniques, making them highly reliable for experiments with a small number of biological replicates.

1. **DESeq2**: Used size factor estimation and a Wald test for significance. Log2 Fold Change (LFC) shrinkage was applied using the `ashr` method to moderate effect sizes for genes with low counts.
2. **edgeR**: Used Trimmed Mean of M-values (TMM) normalization and fitted a negative binomial Generalized Linear Model (GLM). Significance was tested using a Likelihood Ratio Test (LRT).
3. **limma-voom**: Used TMM normalization factors, `voom` transformation to compute observational weights, and fitted a linear model using `lmFit` followed by `eBayes`.

## Thresholds for Significance
To determine significantly differentially expressed genes in each pipeline, the following thresholds were applied:
* Adjusted p-value (FDR) < 0.05
* Absolute Log2 Fold Change ($\|LFC\|$) $\geq$ 1.0 (indicating at least a 2-fold change)

## Key Findings & Consensus
We evaluated the overlap of significant DEGs across all three pipelines. The analysis yielded:
* **DESeq2**: 33 significant genes
* **edgeR**: 8 significant genes
* **limma-voom**: 0 significant genes

Because limma-voom returned 0 genes passing the strict thresholds (a known conservative behavior for datasets with small sample sizes and high biological variance), the strict 3-way consensus is 0.

### Best Pair of Tools
We can logically select the best pair of tools by evaluating their practical utility and agreement on this specific dataset:
1. **limma-voom** is highly conservative on this small-sample dataset (N=3 paired) and returned 0 significant DEGs. Because it fails to nominate any candidates under these strict thresholds, it provides no actionable targets for downstream analysis and is therefore excluded.
2. **DESeq2** nominated 33 DEGs, and **edgeR** nominated 8 DEGs.
3. All 8 DEGs identified by edgeR are completely contained within the 33 DEGs identified by DESeq2. 

Because DESeq2 and edgeR perfectly validate each other on those 8 targets, they form the most robust and practical pair of tools for this task. Focusing on their 2-way intersection, **8 consensus DEGs** were identified:

| Gene Name | Ensembl ID | Chromosome | Biotype | Direction in SARS |
| :--- | :--- | :--- | :--- | :--- |
| **OAS2** | ENSG00000111335 | 12 | protein_coding | Upregulated |
| **IFI27** | ENSG00000165949 | 14 | protein_coding | Upregulated |
| **GPRASP1** | ENSG00000198932 | X | protein_coding | Upregulated |
| **SERPINB2** | ENSG00000197632 | 18 | protein_coding | Downregulated |
| **CXCL11** | ENSG00000169248 | 4 | protein_coding | Upregulated |
| **AC026801.2** | ENSG00000272323 | 5 | lncRNA | Upregulated |
| **LPL** | ENSG00000175445 | 8 | protein_coding | Upregulated |
| **XAF1** | ENSG00000132530 | 17 | protein_coding | Upregulated |

## Chromosome 1 and 2 Candidates
A subset of priority target candidates on Chromosomes 1 and 2 was selected for potential downstream evaluation (Tasks 2 and 3). 
These were filtered for `protein_coding` biotype and relaxed significance thresholds (padj < 0.10 in either pipeline). CRISPR modalities were assigned based on the direction of differential expression (e.g., *CRISPRi (Repression)* for upregulated targets and *CRISPRa (Activation)* for downregulated targets). Notable genes in this list include **GBP4** (Chr 1), **REG3A** (Chr 2), and **SLC16A1** (Chr 1).
