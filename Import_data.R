# Import Data
# 6/29/26
# E. Lamont

# ####################################################### #
##################### LOAD PACKAGES #######################

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
library(scales) # For comma()

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

# ####################################################### #
############### IMPORT PIPELINE SUMMARY DATA ##############

# Marm_1_LRF
Marm_1_pipeSummary <- read.csv("Data/Marm_1_LRF/Pipeline.Summary.Details.csv") %>% 
  select(-X) %>%
  mutate(Run = "Marm_1")

# Marm_2_LRF
Marm_2_pipeSummary <- read.csv("Data/Marm_2_LRF/Pipeline.Summary.Details.csv") %>% 
  select(-X) %>%
  mutate(Run = "Marm_2")

# Marm_3_LRF
Marm_3_pipeSummary <- read.csv("Data/Marm_3_LRF/Pipeline.Summary.Details.csv") %>% 
  select(-X) %>%
  mutate(Run = "Marm_3")

# Marm_4_LRF
Marm_4_pipeSummary <- read.csv("Data/Marm_4_LRF/Pipeline.Summary.Details.csv") %>% 
  select(-X) %>%
  mutate(Run = "Marm_4")

# Ella's H37Rv_Log samples
H37Rv_pipeSummary <- read.csv("Data/PredictTB_Run6/Pipeline.Summary.Details.csv") %>% 
  select(-Outcome, -Arm, -Patient, -X, -Week) %>% 
  filter(Type2 == "H37Rv Log") %>%
  mutate(Run = "PredictTB_Run6")

# Merge the pipeSummaries
All_pipeSummary <- merge(Marm_1_pipeSummary, Marm_2_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, Marm_3_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, Marm_4_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, H37Rv_pipeSummary, all = T)


# Make a second SampleID column
# Remove _S* From names
All_pipeSummary$SampleID2 <- gsub(x = All_pipeSummary$SampleID, pattern = "_S.*", replacement = "")

# Add run numbers to the H37Rvs
All_pipeSummary <- All_pipeSummary %>%
  mutate(SampleID2 = if_else(SampleID2 %in% c("H37Rv_1", "H37Rv_2", "H37Rv_3") 
                             & str_detect(Run, "^Marm_[1234]$"), 
                             paste0(SampleID2, "_Run", str_extract(Run, "[1234]")),
                             SampleID2)) 

# Remove the undetermined
All_pipeSummary <- All_pipeSummary %>%
  filter(SampleID != "Undetermined_S0")

# Add the Binary outcome and rows for extra samples
Binary_Outcome <- read.csv("Data/feature_selection_meta.csv") %>%
  mutate(SampleID2 = base::sub("_[^_]*$", "", SampleID2))
All_pipeSummary <- All_pipeSummary %>%
  full_join(Binary_Outcome %>% select(SampleID2, Outcome2),
            by = "SampleID2")


# ####################################################### #
################## IMPORT SAMPLE METADATA #################

# I edited this in excel to make it better for import
All_metadata <- read.csv("Data/Sample_Metadata/Marm_seq_meta.csv")

All_pipeSummary <- All_pipeSummary %>%
  left_join(All_metadata, by = join_by(Run, SampleID2))

# Add column with JM vs LRF samples
# All_pipeSummary <- All_pipeSummary %>% 
#   mutate(Handler = if_else(str_detect(SampleID2, "JM"), "JM", "LRF")) %>%
#   mutate(Handler = if_else(str_detect(SampleID2, "Rv_Log"), "EIL", Handler))
All_pipeSummary <- All_pipeSummary %>%
  mutate(Handler = case_when(str_detect(SampleID2, "Rv_Log") ~ "EIL",
                             str_detect(SampleID2, "JM") ~ "JM",
                             str_detect(SampleID2, "B") ~ "JM",
                             TRUE ~ "LRF"))

# ####################################################### #
############# ADD MORE COLUMNS TO PIPE SUMMARY ############

All_pipeSummary <- All_pipeSummary %>% 
  mutate(Type = if_else(str_detect(SampleID2, "Marm"), "Marmoset", "H37Rv")) 

All_pipeSummary <- All_pipeSummary %>%
  mutate(across(where(is.character), ~ if_else(str_detect(SampleID2, "Rv_Log") & is.na(.x), "H37Rv", .x)))




# ####################################################### #
############### IMPORT AND PROCESS RAW READS ##############

Run1_RawReads <- read.csv("Data/Marm_1_LRF/Mtb.Expression.Gene.Data.readsM.csv")
Run1_RawReads <- Run1_RawReads %>% 
  dplyr::select(-Undetermined_S0) %>%
  dplyr::rename_with(function(x) gsub("(_S[0-9]+)$", "_Run1\\1", x), 
                     .cols = matches("^H37Rv_[123]_S[0-9]+$"))

Run2_RawReads <- read.csv("Data/Marm_2_LRF/Mtb.Expression.Gene.Data.readsM.csv")
Run2_RawReads <- Run2_RawReads %>% 
  dplyr::select(-Undetermined_S0) %>%
  dplyr::rename_with(function(x) gsub("(_S[0-9]+)$", "_Run2\\1", x), 
                     .cols = matches("^H37Rv_[123]_S[0-9]+$"))

Run3_RawReads <- read.csv("Data/Marm_3_LRF/Mtb.Expression.Gene.Data.readsM.csv")
Run3_RawReads <- Run3_RawReads %>% 
  dplyr::select(-Undetermined_S0) %>%
  dplyr::rename_with(function(x) gsub("(_S[0-9]+)$", "_Run3\\1", x), 
                     .cols = matches("^H37Rv_[123]_S[0-9]+$"))

Run4_RawReads <- read.csv("Data/Marm_4_LRF/Mtb.Expression.Gene.Data.readsM.csv")
Run4_RawReads <- Run4_RawReads %>% 
  dplyr::select(-Undetermined_S0) %>%
  dplyr::rename_with(function(x) gsub("(_S[0-9]+)$", "_Run4\\1", x), 
                     .cols = matches("^H37Rv_[123]_S[0-9]+$"))

JM_RawReads <- read.csv("Data/counts_matrix_JM.csv")
JM_RawReads <- JM_RawReads %>% 
  dplyr::rename_with(function(x) gsub("(_S[0-9]+)$", "_Run4\\1", x), 
                     .cols = matches("^H37Rv_[123]_S[0-9]+$"))

# Ella's H37Rv_Log samples
H37Rv_RawReads <- read.csv("Data/PredictTB_Run6/Mtb.Expression.Gene.Data.readsM.csv") %>% 
  dplyr::select(X, contains("Rv_Log"))

# Merge the RawReads
All_RawReads <- merge(Run1_RawReads, Run2_RawReads, all = T)
All_RawReads <- merge(All_RawReads, Run3_RawReads, all = T)
All_RawReads <- merge(All_RawReads, Run4_RawReads, all = T)
All_RawReads <- merge(All_RawReads, JM_RawReads, all = T)
All_RawReads <- merge(All_RawReads, H37Rv_RawReads, all = T)


# Remove the _S at the end
names(All_RawReads) = gsub(pattern = "_S[0-9]+$", replacement = "", x = names(All_RawReads))

# Keep only the protein coding Rv genes
All_RawReadsf <- All_RawReads %>%
  filter(grepl("^Rv[0-9]+[A-Za-z]?$", X))

All_RawReadsf <- All_RawReadsf %>%
  column_to_rownames("X")

# ####################################################### #
######## CALCULATE TXN COVERAGE FROM Rv GENES ONLY ########

# Count, for each column (sample), how many genes have >= 10 reads
NumGoodReads <- colSums(All_RawReadsf >= 10)

# Add as new column in All_pipeSummary, matching by SampleID2
All_pipeSummary$AtLeast.10.Readsf <- NumGoodReads[All_pipeSummary$SampleID2]

# Add transcriptional coverage
All_pipeSummary <- All_pipeSummary %>% mutate(Txn_Coveragef = round(AtLeast.10.Readsf/4030*100))

# ####################################################### #
#################### VSTB NORMALIZATION ###################
# Blinded VST (VSTB) for PCA. 

# # Keep Genes with at least 10 reads total across all samples
# keep <- rowSums(All_RawReads_f) >= 10
# All_RawReads_f2 <- All_RawReads_f[keep,] # Now only 4025 genes
# 
# # Generate a matrix of integers
# All_RawReads_f2_int <- All_RawReads_f2 %>%
#   mutate(across(everything(), round)) %>% 
#   as.matrix()
# 
# # Normalize without metadata (Blinded)
# All_VSTB <- varianceStabilizingTransformation(All_RawReads_f2_int, fitType = "parametric")
# All_VSTB <- as.data.frame(All_VSTB)

###########################################################
############### COMBATSEQ BATCH CORRECTION ################
# Using CombatSeq

# Keep Genes with at least 10 reads total across all samples (Already run above)
# keep <- rowSums(All_RawReadsf) >= 10
# All_RawReadsf2 <- All_RawReadsf[keep,] # Now only 4025 genes

# # Generate a matrix of integers
# All_RawReadsf2_int <- All_RawReads_f2 %>%
#   mutate(across(everything(), round)) %>%
#   as.matrix()
# 
# count_matrix <- as.matrix(All_RawReadsf2_int)
# # Ensure integer counts
# mode(count_matrix) <- "integer"
# 
# # Reorder metadata to match column order in count_matrix
# meta <- All_pipeSummary[match(colnames(count_matrix), All_pipeSummary$SampleID2), ]
# 
# # Check alignment
# all(meta$SampleID2 == colnames(count_matrix))  # should be TRUE
# # # IF NOT TRUE RUN THESE:
# # ## This will show the samples in count_matrix that don't match metadata
# # # colnames(count_matrix)[!colnames(count_matrix) %in% All_pipeSummary$SampleID2]
# # ## And the opposite: metadata samples not in count_matrix
# # # All_pipeSummary$SampleID2[!All_pipeSummary$SampleID2 %in% colnames(count_matrix)]
# # 
# # 
# # Extract batch (run) and condition (for checking later)
# batch <- meta$Run
# condition <- meta$Type
# # 
# # Run ComBat-Seq
# combat_counts <- ComBat_seq(
#   count_matrix,
#   batch = batch,
#   group = condition # optional, helps preserve biological signal
# )


###########################################################
################ VSTB FROM BATCH CORRECTED ################

# # Normalize without metadata (Blinded)
# All_BC_VSTB <- varianceStabilizingTransformation(combat_counts, fitType = "parametric")
# All_BC_VSTB <- as.data.frame(All_BC_VSTB)


###########################################################
################ MAKE TPM FROM Rv RAW READS ###############

# Bob's TPM includes all the non-coding RNAs, make a new TPM from just the protein-coding genes

source("Function_CalculateTPM.R")
All_tpmf <- CalculateTPM_RvOnly(All_RawReadsf %>% rownames_to_column("X"))

# Remove the _S at the end
names(All_tpmf) = gsub(pattern = "_S[0-9]+$", replacement = "", x = names(All_tpmf))

All_tpmf_log2 <- All_tpmf %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log2(.x))) # Log transform the values

# ####################################################### #
#################### REMOVE DUPLICATES ####################
# Removing duplicates that would make it above the threshold (leaving in the lower ones because they are filtered out later)

my_pipeSummary <- All_pipeSummary %>%
  filter(!SampleID2 %in% c("Marm_LRF_2_re", "Marm_LRF_25_re", "Marm_LRF_42", "Marm_LRF_56_re"))


# ####################################################### #
################### FILTER GOODSAMPLES60 ##################

GoodSamples60_pipeSummary <- my_pipeSummary %>%
  filter(Txn_Coveragef >= 60) %>%
  filter(N_Genomic >= 700000)

GoodSampleList60 <- GoodSamples60_pipeSummary %>%  
  pull(SampleID2) # 71 samples

GoodSamples60_RawReadsf <- All_RawReadsf %>% 
  dplyr::select(all_of(GoodSampleList60))

GoodSamples60_log2tpmf <- All_RawReadsf %>% 
  dplyr::select(all_of(GoodSampleList60))

# GoodSamples60_VSTB <- All_VSTB %>% 
#   dplyr::select(all_of(GoodSampleList60))



###########################################################
##################### SUMMARY NUMBERS #####################

GoodSamples60_pipeSummary %>%
  filter(N_Genomic >= 700000) %>%
  filter(Txn_Coveragef >=60) %>%
  group_by(Type, Cavity_score) %>%
  summarize(N_samples = n())

All_pipeSummary %>%
  filter(N_Genomic < 700000 | Txn_Coveragef <60) %>%
  # filter(Txn_Coverage_f < 60) %>%
  group_by(Type, Cavity_score) %>%
  summarize(N_samples = n())

# ####################################################### #
################### CLEAN UP ENVIRONMENT ##################

# rm(list = ls(pattern = "^tmp"))

# Remove original pipeSummaries
rm(Marm_1_pipeSummary, Marm_2_pipeSummary, Marm_3_pipeSummary, Marm_4_pipeSummary, H37Rv_pipeSummary)

# Remove original raw reads
rm(Run1_RawReads, Run2_RawReads, Run3_RawReads, Run4_RawReads, H37Rv_RawReads)

# Remove VST intermediates
rm(All_RawReads_f2, All_RawReads_f2_int, keep)

# remove batch correction intermediates
# rm(combat_counts, count_matrix, meta, batch, condition)

# Remove random other things
rm(NumGoodReads)


