# PCA
# E. Lamont
# 6/29/26

source("Import_data.R")

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=12),
        # legend.title = element_text(size = 14),
        legend.title = element_text(size=14),
        plot.title = element_text(size=12), 
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(angle = 0, size=12, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12), 
        plot.subtitle = element_text(size=9))

Run_colors <- c(`Marm_1` = "#49006A", 
                `Marm_2` = "#AE017E", 
                `Marm_3`= "#F768A1",
                `Marm_4` = "#FA9FB5",
                `PredictTB_Run6` = "#FCC5C0")

CavityScore_colors <- c("H37Rv" = "#999999",
                        "Normal" = "#A8D5BA",
                        "0: necrotic"  = "#FEE5D9",
                        "0: fibrotic" = "#FC9272", 
                        "1" = "#F4A6A6",
                        "2" = "#EF6C6C",
                        "3" = "#C62828",
                        "5" = "#8E1B1B")

Binary_colors <- c(`H37Rv` = "#999999", 
                   `1` = "#D55E00", 
                   `0`= "#0072B2")

Days_colors <- c(
  "H37Rv"     = "#999999",
  "untreated" = "#E7298A",
  "92"  = "#A6CEE3",
  "136" = "#1F78B4",
  "156" = "#08306B",
  "130" = "#B2DF8A",
  "87"  = "#33A02C",
  "168" = "#006D2C",
  "66"  = "#FDBF6F",
  "73"  = "#FF7F00",
  "48"  = "#E31A1C",
  "72"  = "#FB9A99",
  "84"  = "#6A3D9A")

Tissue_colors <- c(
  "H37Rv"   = "#999999",
  "RUL"     = "#1D91C0", 
  "RML"     = "#41B6C4",
  "RLL"     = "#0C2C84",
  "LUL"     = "#ADDD8E",  
  "LML"     = "#41AB5D",
  "LLL"     = "#005A32",
  "Liver"   = "darkred",
  "ACC (?)" = "#FDBF6F")

Binary_colors2 <- c("#D55E00", "#0072B2")

Type_colors <- c("Marmoset" = "#AE017E",
                 "H37Rv" = "#999999")

###########################################################
#################### ALL SAMPLES VSTB #####################

# Convert gene column to rownames
tmp_data <- All_VSTB 

# Transform the data
tmp_data_t <- as.data.frame(t(tmp_data)) # or my_tpm2

# Make the actual PCA
tmp_PCA <- prcomp(tmp_data_t, scale = F) # Scale is F here because the data is already normalized

# See the % Variance explained
summary(tmp_PCA)
tmp_summary_PCA <- format(round(as.data.frame(summary(tmp_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
tmp_summary_PCA[1,1] # PC1 explains 13.9% of variance
tmp_summary_PCA[2,1] # PC2 explains 4.9% of variance
tmp_summary_PCA[3,1] # PC3 explains 3.3% of variance

# MAKE PCA PLOT with GGPLOT 
tmp_PCA_df <- as.data.frame(tmp_PCA$x[, 1:3]) # Extract the first 3 PCs
tmp_PCA_df <- data.frame(SampleID2 = row.names(tmp_PCA_df), tmp_PCA_df)
tmp_PCA_df <- merge(tmp_PCA_df, All_pipeSummary, by = "SampleID2")


PCA_fig1 <- tmp_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Run)) + 
  geom_point(size = 3.5, alpha = 0.8, stroke = 0.8, shape = 21) +
  geom_text_repel(aes(label = SampleID2), size = 2) + 
  scale_fill_manual(values = Run_colors) +  
  # scale_shape_manual(values = my_fav_shapes) + 
  # geom_text_repel(aes(label = comma(N_Genomic)), size= 2, box.padding = 0.4, segment.color = "grey42", max.overlaps = Inf) + 
  labs(title = "Marm Runs1-4 VSTB",
       subtitle = "All samples",
       x = paste0("PC1: ", tmp_summary_PCA[1,1], "%"),
       y = paste0("PC2: ", tmp_summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_fig1
# ggsave(PCA_fig1,
#        file = paste0("All_VSTB_v1.pdf"),
#        path = "Figures/PCA",
#        # dpi = 150,
#        width = 9, height = 6, units = "in")

PCA_fig2 <- tmp_PCA_df %>%
  ggplot(aes(x = PC2, y = PC3, fill = Run)) +
  geom_point(size = 3.5, alpha = 0.8, stroke = 0.8, shape = 21) +
  scale_fill_manual(values = Run_colors) +
  labs(title = "Marm Runs1-4 VSTB",
       subtitle = "All samples",
       x = paste0("PC2: ", tmp_summary_PCA[2,1], "%"),
       y = paste0("PC3: ", tmp_summary_PCA[3,1], "%")) +
  my_plot_themes
PCA_fig2
# ggsave(PCA_fig2,
#        file = paste0("All_VSTB_v2.pdf"),
#        path = "Figures/PCA",
#        # dpi = 150,
#        width = 9, height = 6, units = "in")

# Batch effect when seeing PC2 vs PC3?


# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Run,
                  colors = my_run_colors)
PCA_3D


# ###################################################### #
################## GOODSAMPLES60 VSTB ####################

# Convert gene column to rownames
tmp_data <- GoodSamples60_VSTB 

# Transform the data
tmp_data_t <- as.data.frame(t(tmp_data)) # or my_tpm2

# Make the actual PCA
tmp_PCA <- prcomp(tmp_data_t, scale = F) # Scale is F here because the data is already normalized

# See the % Variance explained
summary(tmp_PCA)
tmp_summary_PCA <- format(round(as.data.frame(summary(tmp_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
tmp_summary_PCA[1,1] # PC1 explains 18.1% of variance
tmp_summary_PCA[2,1] # PC2 explains 5.1% of variance
tmp_summary_PCA[3,1] # PC3 explains 3.3% of variance

# MAKE PCA PLOT with GGPLOT 
tmp_PCA_df <- as.data.frame(tmp_PCA$x[, 1:3]) # Extract the first 3 PCs
tmp_PCA_df <- data.frame(SampleID2 = row.names(tmp_PCA_df), tmp_PCA_df)
tmp_PCA_df <- merge(tmp_PCA_df, GoodSamples60_pipeSummary, by = "SampleID2")


PCA_fig1 <- tmp_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Run)) + 
  geom_point(size = 3.5, alpha = 0.8, stroke = 0.8, shape = 21) +
  geom_text_repel(aes(label = SampleID2), size = 2) + 
  scale_fill_manual(values = Run_colors) +  
  # scale_shape_manual(values = my_fav_shapes) + 
  # geom_text_repel(aes(label = comma(N_Genomic)), size= 2, box.padding = 0.4, segment.color = "grey42", max.overlaps = Inf) + 
  labs(title = "Marm Runs1-4 VSTB",
       subtitle = "GoodSamples60_VSTB",
       x = paste0("PC1: ", tmp_summary_PCA[1,1], "%"),
       y = paste0("PC2: ", tmp_summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_fig1
ggsave(PCA_fig1,
       file = paste0("GoodSamples60_VSTB_Run_v1.pdf"),
       path = "Figures/PCA",
       # dpi = 600,
       width = 9, height = 6, units = "in")

PCA_fig2 <- tmp_PCA_df %>% 
  ggplot(aes(x = PC2, y = PC3, fill = Run)) + 
  geom_point(size = 3.5, alpha = 0.8, stroke = 0.8, shape = 21) +
  geom_text_repel(aes(label = SampleID2), size = 2) + 
  scale_fill_manual(values = Run_colors) +  
  labs(title = "Marm Runs1-4 VSTB",
       subtitle = "GoodSamples60_VSTB",
       x = paste0("PC2: ", tmp_summary_PCA[2,1], "%"),
       y = paste0("PC3: ", tmp_summary_PCA[3,1], "%")) +
  my_plot_themes
PCA_fig2
ggsave(PCA_fig2,
       file = paste0("GoodSamples60_VSTB_Run_v2.pdf"),
       path = "Figures/PCA",
       # dpi = 600,
       width = 9, height = 6, units = "in")

# Batch effect when seeing PC2 vs PC3?


# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Run,
                  colors = my_run_colors)
PCA_3D

PCA_fig3 <- tmp_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type)) + 
  geom_point(size = 3.5, alpha = 0.8, stroke = 0.8, shape = 21) +
  # geom_text_repel(aes(label = Run), size = 2) + 
  scale_fill_manual(values = Type_colors) +  
  labs(title = "Marm Runs1-4 VSTB",
       subtitle = "GoodSamples60_VSTB",
       x = paste0("PC1: ", tmp_summary_PCA[1,1], "%"),
       y = paste0("PC2: ", tmp_summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_fig3
ggsave(PCA_fig3,
       file = paste0("GoodSamples60_VSTB_Type_v1.pdf"),
       path = "Figures/PCA",
       # dpi = 600,
       width = 9, height = 6, units = "in")
