# Import Data
# 5/27/26
# E. Lamont

################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(plotly)
library(DESeq2)
# library(ggprism) # for add_pvalue()
# library(rstatix) # for adjust_pvalue
# library(ggpmisc) # https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
library(ggrepel)
# library(pheatmap)
# library(dendextend) # May need this for looking at pheatmap clustering
# library(ggplotify) # To convert pheatmaps to ggplots
# library(corrplot)
# library(ggcorrplot)
library(ggfortify) # To make pca plots with plotly
library(edgeR) # for cpm
library(sva) # For ComBat_seq batch correction
library(stringr)
library(readxl) # To import excel files as dataframes

# DuffyTools
# library(devtools)
# install_github("robertdouglasmorrison/DuffyTools")
# library(DuffyTools)
# install_github("robertdouglasmorrison/DuffyNGS")
# BiocManager::install("robertdouglasmorrison/DuffyTools")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biobase")


# Stop scientific notation
# options(scipen = 999) 
# options(scipen = 0) # To revert back to default

###########################################################
############### IMPORT PIPELINE SUMMARY DATA ##############

# Marm_1_LRF
Marm_1_LRF_pipeSummary <- read.csv("Data/Marm_1_LRF/Pipeline.Summary.Details.csv") %>% 
  select(-X) %>%
  mutate(Run = "Marm_1_LRF")


# Make a second SampleID column
# Remove _S* From names
Marm_1_LRF_pipeSummary$SampleID2 <- gsub(x = Marm_1_LRF_pipeSummary$SampleID, pattern = "_S.*", replacement = "")


###########################################################
############### IMPORT AND PROCESS RAW READS ##############

Run1_RawReads <- read.csv("Data/Marm_1_LRF/Mtb.Expression.Gene.Data.readsM.csv")
Run1_RawReads <- Run1_RawReads %>% 
  dplyr::select(-Undetermined_S0)

# Remove the _S at the end
names(Run1_RawReads) = gsub(pattern = "_S[0-9]+$", replacement = "", x = names(Run1_RawReads))

# Keep only the protein coding Rv genes
Run1_RawReads_f <- Run1_RawReads %>%
  filter(grepl("^Rv[0-9]+[A-Za-z]?$", X))


###########################################################
######## CALCULATE TXN COVERAGE FROM Rv GENES ONLY ########

# Count, for each column (sample), how many genes have >= 10 reads
NumGoodReads <- colSums(Run1_RawReads_f >= 10)

# Add as new column in All_pipeSummary, matching by SampleID2
Marm_1_LRF_pipeSummary$AtLeast.10.Reads_f <- NumGoodReads[Marm_1_LRF_pipeSummary$SampleID2]

# Add transcriptional coverage
Marm_1_LRF_pipeSummary <- Marm_1_LRF_pipeSummary %>% mutate(Txn_Coverage_f = round(AtLeast.10.Reads_f/4030*100))


###########################################################
################### FILTER GOODSAMPLES60 ##################

GoodSamples60_pipeSummary <- Marm_1_LRF_pipeSummary %>%
  filter(Txn_Coverage_f >= 60)

GoodSampleList60 <- GoodSamples60_pipeSummary %>%  
  pull(SampleID2) # 29 samples

GoodSamples60_RawReadsf <- Run1_RawReads_f %>% 
  dplyr::select(X, all_of(GoodSampleList60)) %>% 
  column_to_rownames(var = "X")

###########################################################
#################### VSTB NORMALIZATION ###################
# Blinded VST (VSTB) for PCA. Also I don't have the metadata


# Keep genes with >5 counts in at least 50% of samples 
# Not doing this because it threw a parametric error when VST normalizing
# keep <- rowSums(Run1_RawReads_f > 5) >= 0.5 * ncol(Run1_RawReads_f)
# Run1_RawReadsf2 <- Run1_RawReads_f[keep, ] # now only 3707 genes (instead of 4030)

# Keep Genes with at least 10 reads total across all samples
keep <- rowSums(GoodSamples60_RawReadsf) >= 10
GoodSamples60_RawReadsf2 <- GoodSamples60_RawReadsf[keep,] # Now only 4025 genes

# Generate a matrix of integers
GoodSamples60_RawReadsf2_int <- GoodSamples60_RawReadsf2 %>%
  mutate(across(everything(), round)) %>% 
  as.matrix()

# Normalize without metadata (Blinded)
GoodSamples60_VSTB <- varianceStabilizingTransformation(GoodSamples60_RawReadsf2_int, fitType = "parametric")
GoodSamples60_VSTB <- as.data.frame(GoodSamples60_VSTB)




