# PCA
# E. Lamont
# 5/27/26

source("Import_data.R")

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9))

###########################################################
##################### PCA RUN 1 VSTB ######################

# Convert gene column to rownames
my_data <- GoodSamples60_VSTB 

# Transform the data
my_data_t <- as.data.frame(t(my_data)) # or my_tpm2

# Make the actual PCA
my_PCA <- prcomp(my_data_t, scale = F) # Scale is F here because the data is already normalized

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 12.1% of variance
summary_PCA[2,1] # PC2 explains 8.9% of variance
summary_PCA[3,1] # PC3 explains 7.2% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
# my_PCA_df <- merge(my_PCA_df, GoodSamples60_pipeSummary, by = "SampleID2")

PCA_fig <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(size = 5, alpha = 0.8, stroke = 0.8) +
  # geom_text_repel(aes(label = SampleID), size = 2.5) + 
  # scale_fill_manual(values = my_fav_colors) +  
  # scale_shape_manual(values = my_fav_shapes) + 
  geom_text_repel(aes(label = SampleID), size= 2, box.padding = 0.4, segment.color = "black", max.overlaps = Inf) + 
  labs(title = "Marm_1 VST Blinded",
       subtitle = "GoodSamples60_RawReadsf2 -> DESeq2 VST blinded",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_fig

# PCA_fig <- my_PCA_df %>% 
#   ggplot(aes(x = PC1, y = PC2, fill = Type2, shape = Type2)) + 
#   geom_point(aes(fill = Type2, shape = Type2), size = 5, alpha = 0.8, stroke = 0.8) +
#   # geom_text_repel(aes(label = Lineage), size = 2.5) + 
#   scale_fill_manual(values = my_fav_colors) +  
#   scale_shape_manual(values = my_fav_shapes) + 
#   # geom_text_repel(aes(label = Patient), size= 2, box.padding = 0.4, segment.color = "black", max.overlaps = Inf) + 
#   labs(title = "GoodSputum60 (Run1-4) VST Blinded",
#        subtitle = "RawReadf -> DESeq2 VST blinded",
#        x = paste0("PC1: ", summary_PCA[1,1], "%"),
#        y = paste0("PC2: ", summary_PCA[2,1], "%")) +
#   my_plot_themes
# PCA_fig

