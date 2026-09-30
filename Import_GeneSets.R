# Import Gene Sets that have been made elsewhere
# E. Lamont
# 7/29/25


###########################################################
###################### LOAD GENE SETS #####################

# To put them in a list of lists
rda_files <- list.files("Data/GeneSet_Data", pattern = "\\.rda$", full.names = TRUE)
allGeneSetList <- list()

for(file in rda_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  env <- new.env()
  load(file, envir = env)  # loads allGeneSets into env
  allGeneSetList[[file_name]] <- env$allGeneSets  # store it in our list
  rm(env)
}

# Now update the names for each gene set in the list of lists:
allGeneSetList <- lapply(allGeneSetList, function(gset) {
  if (!is.null(gset)) {
    names(gset) <- gsub("<.*", "", names(gset))
  }
  return(gset)
})



###########################################################
################# TO MAKE NEW GENE SETS ###################

# First! Make the .csv file by hand. Needs a "Gene" and a "GeneSet" column

# XXX_GeneSets <- read.csv("Data/GeneSet_Data/XXX_GeneSets.csv")

# Create a list where each GeneSet is a named element
# Needs to be called allGeneSets so it is easier to load with all the others
# allGeneSets <- split(XXX_GeneSets$Gene, XXX_GeneSets$GeneSet)

# SAVE AS RDA FOR LATER 
# save(allGeneSets, file = "Data/GeneSet_Data/XXXGeneSets.rda")

# 8/24/26: Lance Domains
# LanceDomains <- read.csv("Data/GeneSet_Data/LanceDomains_Aug2026.csv")
# allGeneSets <- split(LanceDomains$Gene, LanceDomains$GeneSet)
# save(allGeneSets, file = "Data/GeneSet_Data/LanceDomains_Aug2026.rda")

# 8/31/26: Subset gene sets interesting in W2 Relapse vs Cure
# Ella_W2Subsets_2026.08.31 <- read.csv("Data/GeneSet_Data/Ella_W2Subsets_2026.08.31.csv")
# allGeneSets <- split(Ella_W2Subsets_2026.08.31$Gene, Ella_W2Subsets_2026.08.31$GeneSet)
# save(allGeneSets, file = "Data/GeneSet_Data/Ella_W2Subsets_2026.08.31.rda")

###########################################################
############ iMODULONS: MAKE LISTS OF GROUPS ##############

# CentralCarbon_iModulons <- c("Peptidoglycan Biosynthesis", "Central Carbon Metabolism", "Fumarate Reductase", "PrpR", "BkaR", "Nicotinate Metabolism") # Needs the \\ to detect it
# CentralCarbon_iModulons_pattern <- str_c(CentralCarbon_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or
# 
# AminoAcid_iModulons <- c("GroEL-GroES Complex", "Leucine Related", "LysG", "ArgR") # Needs the \\ to detect it
# AminoAcid_iModulons_pattern <- str_c(AminoAcid_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or
# 
# NucleicAcid_iModulons <- c("PyrR", "Rv0135\\+Rv1019", "Nucleic Acid Hydrolysis") # Needs the \\ to detect it
# NucleicAcid_iModulons_pattern <- str_c(NucleicAcid_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or
# 
# FattyAcid.Cholesterol_iModulons <- c("Fatty Acid Biosynthesis", "KstR2", "Mycofactocin Synthesis Pathway", "FasR", "Polyketide Synthase Complex", "Rv0681") # Needs the \\ to detect it
# FattyAcid.Cholesterol_iModulons_pattern <- str_c(FattyAcid.Cholesterol_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or
# 
# Metal_iModulons <- c("RicR", "IdeR", "M-box", "Zur", "Hpt-2b Induced") # Needs the \\ to detect it
# Metal_iModulons_pattern <- str_c(Metal_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or
# 
# SulfurMetabolism_iModulons <- c("Sulfur Metabolism") # Needs the \\ to detect it
# SulfurMetabolism_iModulons_pattern <- str_c(SulfurMetabolism_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or
# 
# Growth_iModulons <- c("Positive Regulation of Growth") # Needs the \\ to detect it
# Growth_iModulons_pattern <- str_c(Growth_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or
# 
# Redox_iModulons <- c("DevR-1", "WhiB4", "DevR-2", "WhiB1", "WhiB4/IdeR", "Rv1828/SigH", "Rv1776c\\+WhiB4", "VirS", "WhiB6") # Needs the \\ to detect it
# Redox_iModulons_pattern <- str_c(Redox_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or
# 
# AcidStress_iModulons <- c("MarR") # Needs the \\ to detect it
# AcidStress_iModulons_pattern <- str_c(AcidStress_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or
# 
# Antibiotic_iModulons <- c("Lsr2", "Blal", "Rv0078\\+Rv2034", "WhiB7", "IniR") # Needs the \\ to detect it
# Antibiotic_iModulons_pattern <- str_c(Antibiotic_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or
# 
# Virulence.Persistence_iModulons <- c("Rv0576", "Mce1R", "SigH", "PhoP", "Mce3R", "MprA", "PDIM\\;PGL Synthesis", "Rv2488c", "SigC", "SigD", "MbcA\\+Rv3249c\\+Rv3066", "SigK") # Needs the \\ to detect it
# Virulence.Persistence_iModulons_pattern <- str_c(Virulence.Persistence_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or


###########################################################
################### DOWNLOAD GENE SETS ####################

# allGeneSetList$MTb.iModulons

# iModulons_df <- as.data.frame(allGeneSetList$MTb.iModulons)
# write.csv(iModulons_df, file = "iModulons_GeneSets.csv")

# iModulons_df <- data.frame(
#   iModulon = names(allGeneSetList$MTb.iModulons),
#   genes = sapply(allGeneSetList$MTb.iModulons, paste, collapse = ", ")
# )
# # 
# write.csv(iModulons_df, file = "Data/iModulons_GeneSets_2.csv", row.names = FALSE)

# iModulons_df2 <- iModulons_df %>% separate_longer_delim(genes, delim = ", ") %>%
#   rename(GeneSet = iModulon, Gene = genes) %>%
#   mutate(Reference = "iModulons") %>% 
#   mutate(iModulonCategory = case_when(
#     str_detect(GeneSet, CentralCarbon_iModulons_pattern) ~ "Central Carbon",
#     str_detect(GeneSet, AminoAcid_iModulons_pattern) ~ "Amino Acid",
#     str_detect(GeneSet, NucleicAcid_iModulons_pattern) ~ "Nucleic Acid",
#     str_detect(GeneSet, FattyAcid.Cholesterol_iModulons_pattern) ~ "Fatty Acid_Cholesterol",
#     str_detect(GeneSet, Metal_iModulons_pattern) ~ "Metal",
#     str_detect(GeneSet, SulfurMetabolism_iModulons_pattern) ~ "Sulfur",
#     str_detect(GeneSet, Growth_iModulons_pattern) ~ "Growth",
#     str_detect(GeneSet, Redox_iModulons_pattern) ~ "Redox",
#     str_detect(GeneSet, AcidStress_iModulons_pattern) ~ "Acid Stress",
#     str_detect(GeneSet, Antibiotic_iModulons_pattern) ~ "Antibiotic",
#     str_detect(GeneSet, Virulence.Persistence_iModulons_pattern) ~ "Virulence_Persistence",
#     TRUE ~ "Other"
#   ))

# write.csv(iModulons_df2, file = "Data/GeneSet_Data/iModulons_GeneSets.csv", row.names = FALSE)


# Walter2015_GeneSets <- data.frame(
#   GeneSet = names(allGeneSetList$Walter2015GeneSets),
#   Gene = sapply(allGeneSetList$Walter2015GeneSets, paste, collapse = ", ")) %>%
#   separate_longer_delim(Gene, delim = ", ") %>%
#   mutate(Reference = "Walter2015")
# write.csv(Walter2015_GeneSets, file = "Data/GeneSet_Data/Walter2015_GeneSets.csv", row.names = FALSE)

# 8/17/26: TFOE Regulons
# TFOE.Regulons_df <- data.frame(
#   GeneSet = names(allGeneSetList$MTb.TFOE.Regulons),
#   Gene = sapply(allGeneSetList$MTb.TFOE.Regulons, paste, collapse = ", ")) %>%
#   separate_rows(Gene, sep = ",\\s*")
# write.csv(TFOE.Regulons_df, file = "Data/GeneSet_Data/TFOE.Regulons.csv", row.names = FALSE)

