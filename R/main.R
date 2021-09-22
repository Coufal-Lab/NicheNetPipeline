# Michael Zakariaie
# Generalized NicheNet Pipeline

library(nichenetr)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(RColorBrewer)

#' Configure function
#' 
#' Reads in a yml file to set working directory, paths to data, and constants such as tpm_cutoff
#' 
#' @param The global path to the config file
#' @return no return value
#' @export
configure <- function(config_file) {
  config <- read_yaml(config_file)
  tpm_cutoff <<- config$tpm_cutoff
  setwd(config$working_directory)
  
  output_directory <<- config$output_directory
  if (!file.exists(output_directory)) dir.create(output_directory)

  reciever_expression <<- read.csv(config$reciever_expression)
  sender_expression <<- read.csv(config$sender_expression)
  
  sender_specific_genes <<- read.csv(config$sender_specific_genes)
  names(sender_specific_genes)[1] <<- 'gene'
}

configure('/Users/michaelzakariaie/Desktop/nichenetConfig.yml')

# These are the ligands/targets you want a figure of the signaling pathway for
ligands_all = c("NLGN3")
targets_all = c("JMJD6")
# ------------------------------------------------------------

# devtools::create('NicheNetPipeline')
# library(NicheNetPipeline)
lr_network <- readRDS(system.file("extdata", "lr_network.rds", package = "NicheNetPipeline"))
ligand_target_matrix <- readRDS(system.file("extdata", "ligand_target_matrix.rds", package = "NicheNetPipeline"))
weighted_networks <- readRDS(system.file("extdata", "weighted_networks.rds", package = "NicheNetPipeline"))

ligand_tf_matrix <- readRDS(system.file("extdata", "ligand_tf_matrix.rds", package = "NicheNetPipeline"))
sig_network <- readRDS(system.file("extdata", "signaling_network.rds", package = "NicheNetPipeline"))
gr_network <- readRDS(system.file("extdata", "gr_network.rds", package = "NicheNetPipeline"))

# Reading in data and setting constants ^^

# Apply TPM cutoff
reciever_expressed_genes = reciever_expression[reciever_expression$fetal_mg_avg > tpm_cutoff, "gene"]
sender_expressed_genes = sender_expression[sender_expression$fetal.avg>tpm_cutoff, "gene"]

# Set geneset_oi
geneset_oi <- sender_specific_genes[sender_specific_genes$gene %in% reciever_expressed_genes, "gene"]
geneset_oi <- geneset_oi[geneset_oi %in% rownames(ligand_target_matrix)]
geneset_oi <- as.vector(geneset_oi)

# Set Background_expressed_genes
background_expressed_genes = reciever_expressed_genes[reciever_expressed_genes %in% rownames(ligand_target_matrix)]

# Get potential ligands
ligands = lr_network %>% pull(from) %>% unique()
expressed_genes_bulk = intersect(ligands,sender_expressed_genes)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,reciever_expressed_genes)

lr_network_expressed = lr_network %>% filter(from %in% expressed_genes_bulk & to %in% expressed_receptors) 

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()


# Run NicheNet
ligand_activities = predict_ligand_activities(geneset = geneset_oi,
                                              background_expressed_genes = background_expressed_genes,
                                              ligand_target_matrix = ligand_target_matrix,
                                              potential_ligands = potential_ligands)

# save ligand activities
ligand_activities %>% arrange(-pearson) %>% write_tsv(path=paste(output_directory, 'ligand_activity.txt', sep="/"))

best_upstream_ligands = ligand_activities %>% top_n(10, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

# examine the overall predictions of ligand activity when deciding how to further analyze data
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange", bins = 30)  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(10, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
png(file=paste(output_directory,"ligand_activity_histogram",sep='/'),
    width=600, height=350)
p_hist_lig_activity
dev.off()

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows()



active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.1)

# ordering to make plotting nice
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized Sender ligands",
                                                                    "Reciever Specific Genes",
                                                                    color = "purple",
                                                                    legend_position = "top",
                                                                    x_axis_position = "top",
                                                                    legend_title = "Regulatory potential") + theme(axis.text.x = element_text(face = "italic"))

png(file=paste(output_directory,"ligand_target_network.png",sep='/'),
    width=1000, height=500)
p_ligand_target_network
dev.off()

#### Examine ligand-receptor network for top ranked ligands ####

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

# plot heatmap of ligand-receptor interactions
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized Cortex-ligands","Receptors expressed by Reciever Cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
png(file=paste(output_directory,"ligand_receptor_network.png",sep='/'),
    width=1000, height=500)
p_ligand_receptor_network
dev.off()


#### Visualize expression of ligands and target genes ####

# prep matrix of ligand activities from prediction step
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")

p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized LSEC-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")
png(file=paste(output_directory,"ligand_pearson_corr.png",sep='/'),
    width=300, height=500)
p_ligand_pearson
dev.off()

figures_without_legend = plot_grid(p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                                   p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                   p_ligand_receptor_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                   p_hist_lig_activity + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                   align = 'hv',
                                   nrow = 2)

legends = plot_grid(
  as_ggplot(get_legend(p_ligand_pearson)),
  as_ggplot(get_legend(p_ligand_target_network)),
  as_ggplot(get_legend(p_ligand_receptor_network)),
  as_ggplot(get_legend(p_hist_lig_activity)),
  nrow = 2,
  align = "h")

all_figures = plot_grid(figures_without_legend,
          legends,
          rel_heights = c(10,2), nrow = 2, align = "hv")
png(file=paste(output_directory,"All_Figures.png",sep='/'),
    width=1500, height=900)
all_figures
dev.off()


# Signaling Path Graphs
active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

signaling_network_graph = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

DiagrammeR::render_graph(signaling_network_graph, layout = "tree")

