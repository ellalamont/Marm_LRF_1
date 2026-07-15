# Use RNASeqPower package to determine power
# E. Lamont 
# 2/3/26



###########################################################
################ COEFFICIENT OF VARIATION #################
# Need this for estimating power

# Coefficient of Variation
# = Standard Deviation / mean

# https://support.bioconductor.org/p/9135351/
# edgeR::estimateCommonDisp

# I guess I do have technical sequencing replicates... maybe I could use them? Not done...

# Try and estimate CV from RawReads 60%Cov
Subset_RawReadsf <- GoodSamples60_RawReadsf %>% 
  # column_to_rownames("X") %>%
  select(any_of(GoodSamples60_pipeSummary$SampleID2)) %>%
  select(-contains("Rv_"))

# Using the EdgeR CV calcutator
y <- DGEList(counts = Subset_RawReadsf)
dge <- estimateCommonDisp(y, verbose = T)
# Disp = 3.34304 , BCV = 1.8284 
sqrt(dge$common.dispersion) # 1.828398

# What if I look at the broth samples
Broth_RawReadsf <- GoodSamples60_RawReadsf %>% 
  select(contains("Rv"))
y3 <- DGEList(counts = Broth_RawReadsf)
dge3 <- estimateCommonDisp(y3, verbose = T)
# Disp = 0.14748 , BCV = 0.384 


###########################################################
################## MEAN SEQUENCING DEPTH ##################

mean_depth <- mean(rowMeans(Subset_RawReadsf))
mean_depth
# 6331.369

###########################################################
####################### RNASeqPower #######################

# https://bioconductor.org/packages/release/bioc/manuals/RNASeqPower/man/RNASeqPower.pdf
# https://github.com/royfrancis/shiny-rnaseq-power

# Shiny web app
# https://rnaseq-power.serve.scilifelab.se/app/rnaseq-power

# The paper
# https://www.liebertpub.com/doi/10.1089/cmb.2012.0283

BiocManager::install("RNASeqPower")
library(RNASeqPower)

# Calculating effect size, which is the detectable fold change in expression

# rnapower(depth = 10, n = 12, n2 = 200, cv = 0.72, alpha = 0.05, power = 0.8)
# [1] 1.924735 should be the effect size (Relative expression effect)

# Best case:
# rnapower(depth = 10, n = 21, n2 = 299, cv = 0.72, alpha = 0.05, power = 0.8)
# [1] 1.64437 effect size

rnapower(depth = 6000, n = 8, cv = 1.8284, alpha = 0.05, power = 0.8)
# [1] 12.9523


# NOT DONE HERE!!!!
# ###########################################################
# ########################### pwr ###########################
# 
# install.packages("pwr")
# library(pwr)
# 
# # This is less appropriate for RNA-seq work! 
# 
# pwr.t2n.test(n1 = 12, n2 = 200, d = NULL, sig.level = 0.05, power = 0.8, alternative = c("two.sided"))
# # t test power calculation 
# # 
# # n1 = 12
# # n2 = 200
# # d = 0.8364906
# # sig.level = 0.05
# # power = 0.8
# # alternative = two.sided
# 
# pwr.t.test(n = NULL, d = 2, sig.level = 0.05, power = 0.8, alternative = c("two.sided"))
# # Two-sample t test power calculation 
# # 
# # n = 5.089995
# # d = 2
# # sig.level = 0.05
# # power = 0.8
# # alternative = two.sided
# # 
# # NOTE: n is number in *each* group
# 
# pwr.t.test(n = NULL, d = 1, sig.level = 0.05, power = 0.8, alternative = c("two.sided"))
# # Two-sample t test power calculation 
# # 
# # n = 16.71472
# # d = 1
# # sig.level = 0.05
# # power = 0.8
# # alternative = two.sided
# # 
# # NOTE: n is number in *each* group
