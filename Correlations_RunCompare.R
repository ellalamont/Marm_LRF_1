# Compare the same sample between different runs

# source("Import_data.R") 


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        # axis.text.x = element_text(angle = 45, size=14, vjust=1, hjust=1),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9))


# ####################################################### # 
############### VSTB H37Rv RUN1 vs RUN2 ###################

# Using all the genes
Sample1 <- "H37Rv_1_Run1" 
Sample2 <- "H37Rv_1_Run4" 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")


# ####################################################### #
######################## VSTB Marm_LRF_10 #################

# Using all the genes
Sample1 <- "Marm_LRF_10" # 
Sample2 <- "Marm_LRF_10_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")

# ####################################################### #
######################## VSTB Marm_LRF_12 #################

# Using all the genes
Sample1 <- "Marm_LRF_12" # 
Sample2 <- "Marm_LRF_12_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")

# ####################################################### #
######################## VSTB Marm_LRF_13 #################

# Using all the genes
Sample1 <- "Marm_LRF_13" # 
Sample2 <- "Marm_LRF_13_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")


# ####################################################### #
######################## VSTB Marm_LRF_14 #################

# Using all the genes
Sample1 <- "Marm_LRF_14" # 
Sample2 <- "Marm_LRF_14_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")


# ####################################################### #
######################## VSTB Marm_LRF_15 #################

# Using all the genes
Sample1 <- "Marm_LRF_15" # 
Sample2 <- "Marm_LRF_15_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")


# ####################################################### #
######################## VSTB Marm_LRF_2 #################

# Using all the genes
Sample1 <- "Marm_LRF_2" # 
Sample2 <- "Marm_LRF_2_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")


# ####################################################### #
######################## VSTB Marm_LRF_20 #################

# Using all the genes
Sample1 <- "Marm_LRF_20" # 
Sample2 <- "Marm_LRF_20_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")


# ####################################################### #
######################## VSTB Marm_LRF_22 #################

# Using all the genes
Sample1 <- "Marm_LRF_22" # 
Sample2 <- "Marm_LRF_22_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")


# ####################################################### #
######################## VSTB Marm_LRF_25 #################

# Using all the genes
Sample1 <- "Marm_LRF_25" # 
Sample2 <- "Marm_LRF_25_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")

Sample1 <- "Marm_LRF_25" # 
Sample2 <- "Marm_LRF_25_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_BC_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "BC_VSTB, Pearson correlation",
       x = paste0("BC_VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("BC_VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_BC_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")

# ####################################################### #
######################## VSTB Marm_LRF_3 #################

# Using all the genes
Sample1 <- "Marm_LRF_3" # 
Sample2 <- "Marm_LRF_3_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")

# ####################################################### #
######################## VSTB Marm_LRF_4 #################

# Using all the genes
Sample1 <- "Marm_LRF_4" # 
Sample2 <- "Marm_LRF_4_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")


# ####################################################### #
######################## VSTB Marm_LRF_42 #################

# Using all the genes
Sample1 <- "Marm_LRF_42" # 
Sample2 <- "Marm_LRF_42_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")


# ####################################################### #
######################## VSTB Marm_LRF_56 #################

# Using all the genes
Sample1 <- "Marm_LRF_56" # 
Sample2 <- "Marm_LRF_56_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")


# ####################################################### #
######################## VSTB Marm_LRF_7 #################

# Using all the genes
Sample1 <- "Marm_LRF_7" # 
Sample2 <- "Marm_LRF_7_re" # 
Sample1_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample1) %>% pull(Txn_Coverage_f)
Sample2_TxnCov <- All_pipeSummary %>% filter(SampleID2 == Sample2) %>% pull(Txn_Coverage_f)
ScatterCorr <- All_VSTB %>% 
  rownames_to_column("Gene") %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "VSTB, Pearson correlation",
       x = paste0("VSTB ", Sample1, " (", Sample1_TxnCov, "% TxnCov)"), 
       y = paste0("VSTB ", Sample2, " (", Sample2_TxnCov, "% TxnCov)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, "_VSTB.pdf"),
       path = "Figures/Correlations_RunCompare",
       width = 7, height = 5, units = "in")
