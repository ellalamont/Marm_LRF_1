# Make graphs of the N_Genomic and P_Genomic
# E. Lamont
# 5/27/26

source("Import_data.R") # All_pipeSummary

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
        # legend.title = element_blank(),
        legend.title = element_text(size=10),
        plot.title = element_text(size=10), 
        # axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        # axis.text.x = element_text(angle = 45, size=10, vjust=1, hjust=1),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=9))


# Stop scientific notation
# options(scipen = 999) 
# options(scipen = 0) # To revert back to default

Run_colors <- c(`Marm_1` = "#49006A", 
                `Marm_2` = "#AE017E", 
                `Marm_3`= "#F768A1")

CavityScore_colors <- c("H37Rv" = "#999999",
                        "Normal" = "#A8D5BA",
                        "0: necrotic"  = "#FEE5D9",
                        "0: fibrotic" = "#FC9272", 
                        "1" = "#F4A6A6",
                        "2" = "#EF6C6C",
                        "3" = "#C62828",
                        "5" = "#8E1B1B")


# "Lesion_ID"          "Animal_ID"           "Cavity_score"        "Only._caseum"        "Only_Cellular"      "Mixed"               "Normal._tissue"      "Vaccinated"          "Days.post.treatment" "Tissue.Location"     "Sibling"             "Strain"              "AtLeast.10.Reads_f" "Txn_Coverage_f"      "Handler"            


# ####################################################### #
################### N_GENOMIC VS TXN COV ##################


fig1 <- All_pipeSummary %>% 
  ggplot(aes(x = N_Genomic, y = Txn_Coverage_f, fill = Cavity_score)) + 
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, shape = 21) +
  geom_text_repel(aes(label = SampleID2), size = 2.5) + 
  scale_fill_manual(values = CavityScore_colors) +  
  scale_x_continuous(trans='log10') +
  geom_hline(yintercept = 60, linetype = "dashed", alpha = 0.5) + 
  geom_vline(xintercept = 700000, linetype = "dashed", alpha = 0.5) + 
  labs(title = "Marm Runs1-3",
       subtitle = "All samples",
       x = paste0("# reads aligning to Mtb"),
       y = paste0("% transcriptional coverage")) +
  my_plot_themes
fig1








