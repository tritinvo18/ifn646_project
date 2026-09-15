

# https://github.com/crisprVerse/Tutorials/tree/master
# https://www.nature.com/articles/s41467-022-34320-7


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# #BiocManager::install(version="release") 
#BiocManager::install("crisprScore")
# 
# if (!requireNamespace("devtools", quietly = TRUE))
#   install.packages("devtools")
# devtools::install_github("crisprVerse/crisprDesignData")
# 
# 
# BiocManager::install("BSgenomeForge")
# BiocManager::install("forgeBSgenomeDataPkgFromNCBI")
# 
# 
# x_genome <- forgeBSgenomeDataPkgFromNCBI()
# 
# BSgenomeForge::forgeBSgenomeDataPkgFromNCBI(assembly_accession="GCF_000001405.26",
#                              pkg_maintainer="info@ncbi.nlm.nih.gov",
#                              organism="Homo sapiens",
#                              destdir=tempdir())
# 
# install.packages("devtools")
# devtools::install_github("crisprVerse/crisprScoreData")

library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Rbowtie)
library(txdbmaker)
library(crisprScore)
library(crisprScoreData)



# Create / load bowtie index
fasta <- r"(C:\Users\hanna\OneDrive\Documents\QUT - Semester 6 202607\IFN646 Biomedical Data Science\Group Project\hg38.fa)"
outdir <- r"(C:\Users\hanna\OneDrive\Documents\QUT - Semester 6 202607\IFN646 Biomedical Data Science\bowtie_index\v2)"
# Rbowtie::bowtie_build(fasta,
#                       outdir=outdir,
#                       force=TRUE,
#                       prefix="hg38")
bowtie_index <- file.path(outdir, "hg38")


data("tss_human", package="crisprDesignData")
data(SpCas9, package="crisprBase")
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

# Build gene annotation
gtf_file_path <- r"(C:\Users\hanna\OneDrive\Documents\QUT - Semester 6 202607\IFN646 Biomedical Data Science\Group Project\GCF_000001405.26_GRCh38_genomic.gtf)"

txdb <- makeTxDbFromGFF(
  gtf_file_path,
  format = "gtf"
)

grList <- TxDb2GRangesList(txdb)
GenomeInfoDb::genome(grList) <- "hg38"


#head(grList)

# Get target region on SERPINB2 gene
target_window <- c(0, 500)
target_region <- queryTss(tss_human,
                          #featureType="cds", #Coding region
                          queryColumn="gene_symbol",
                          queryValue="SERPINB2", #"KRAS", #
                          tss_window=target_window
                          )

target_region

gs <- findSpacers(target_region,
                  crisprNuclease=SpCas9,
                  bsgenome=BSgenome.Hsapiens.UCSC.hg38)


gs <- addSequenceFeatures(gs)

gs <- addSpacerAlignments(gs,
                          txObject=grList,
                          aligner_index=bowtie_index,
                          bsgenome=bsgenome,
                          n_mismatches=1)

#alignments(gs)

gs <- addOffTargetScores(gs)
gs <- addOnTargetScores(gs, methods="crisprater")

# Add search for restriction enzymes
gs <- addRestrictionEnzymes(gs)


# Choose best guides
guideSet <- gs
guideSet <- guideSet[guideSet$percentGC>=20]
guideSet <- guideSet[guideSet$percentGC<=80]
guideSet <- guideSet[!guideSet$polyT]

guideSet

guideSet[guideSet$pam == 'GGG']

guideSet[guideSet$polyA == FALSE & guideSet$polyC == FALSE & guideSet$polyG == FALSE & guideSet$polyT == FALSE]

# Creating an ordering index based on the CRISPRater score:
# Using the negative values to make sure higher scores are ranked first:
o <- order(-guideSet$score_crisprater) 
# Ordering the GuideSet:
guideSet <- guideSet[o]
head(guideSet)

# onTargets(gs, columnName = "alignments", unlist = TRUE, use.names = TRUE)
offTargets(
  gs,
  columnName = "alignments",
  max_mismatches = Inf,
  unlist = TRUE,
  use.names = TRUE
)

