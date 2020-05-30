# Installation of sleuth can be painless in case of using conda package, but R should be added to conda enviroment.

library(sleuth)
library(ggplot2)

# This part will take about 30 Gb of memory. File sample_table_for_sleuth.csv contains information about samples with paths to kallisto results. 
# It is important to have columns "path" and "sample", others are optional and might have other names.

s2c <- read.csv('sample_table_for_sleuth.csv')
s2c$path <- as.character(s2c$path)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, num_cores = 6) # Very resource-consuming step - preperation of sleuth object
so <- sleuth_fit(so, ~gastrula_state, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full') # Several models 
so <- sleuth_wt(so, which_beta = 'gastrula_statebefore_gastrula') # Wald's test require replicates!

# Sleuth object can be saved on remote server, downloaded to local machine and then loaded to script with sleuth_load function
sleuth_save(so, "sleuth_obj") 
so <- sleuth_load("sleuth_obj")

# Downstream step are not resource-consuming and can be done on local machine.

sleuth_live(so) # Nice shiny application for results exploring

models(so) # Available models
tests(so) # Available tests

# Graphs 
pca_graph <- plot_pca(so, color_by = "gastrula_state", size = 5, point_size = 10, pc_x = 2, pc_y = 1, text_labels = FALSE)
ggsave("pca.png", pca_graph)
plot_pc_variance <- plot_pc_variance(so,  units = "tpm")
ggsave("varience.png", plot_pc_variance)
volcano <- plot_volcano(so, test = "gastrula_statebefore_gastrula", test_type = "wald", which_model = "full")
ggsave("fold_change.png", volcano)
loads1 <- plot_loadings(so, units = 'tpm', pc_input = 1, pc_count = 20)
ggsave("loads1.png", loads1)
print(loads1$data)
loads2 <- plot_loadings(so, units = 'tpm', pc_input = 2, pc_count = 20)
ggsave("loads2.png", loads2)
print(loads2$data)

# Creating table with interesting result
sleuth_table <- sleuth_results(so, test = 'gastrula_statebefore_gastrula', test_type = "wt", which_model = "full", show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
sleuth_significant <- order(sleuth_significant, by = 'b') # b column contains information about fold change in expression
most_changed <- rbind(head(sleuth_significant, 25), 
                      tail(sleuth_significant, 25))

# Creating heatmap with interesting results
heatmap <- plot_transcript_heatmap(so, transcripts = most_changed$target_id, cluster_transcripts = TRUE, units = "tpm")
ggsave("heatmap.png", heatmap)



