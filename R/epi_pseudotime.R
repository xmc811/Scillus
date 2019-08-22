
library(monocle)
library(colormap)
library(colorRamps)

# Epithelial pseudotime analysis

epi_cds <- seurat_to_monocle(dataset = gfp_combined, subset = T, clusters = c('Basal 1','Luminal 1','Basal 2','Luminal 2'))

epi_cds <- analyze_monocle(epi_cds, rev = F)

#Visualization

svglite(file = "./figures/epi_trajectory_group.svg", width = 10, height = 10)
plot_cell_trajectory(epi_cds, color_by = "group", show_branch_points = F) + scale_color_manual(values = c('#7bccc4','#f03b20'))
dev.off()

svglite(file = "./figures/epi_trajectory_group_split.svg", width = 5, height = 10)
plot_cell_trajectory(epi_cds, color_by = "group", show_branch_points = F) + scale_color_manual(values = c('#7bccc4','#f03b20')) +
        facet_wrap(~group, nrow = 2)
dev.off()

svglite(file = "./figures/epi_trajectory_cluster.svg", width = 10, height = 10)
plot_cell_trajectory(epi_cds, color_by = "cluster", show_branch_points = F) + 
        scale_color_manual(name = '', values = c('#f03b20','#a1d99b','#6baed6','#fed976'))
dev.off()

svglite(file = "./figures/epi_trajectory_cluster_split.svg", width = 10, height = 7)
plot_cell_trajectory(epi_cds, color_by = "group", show_branch_points = F) + 
        scale_color_manual(values = c('#a1d99b','#f03b20'), 
                           labels = c('Early','Advanced'),
                           name = 'Stage') +
        theme(title = element_blank()) +
        facet_wrap(~cluster + group, nrow = 2)
dev.off()

svglite(file = "./figures/epi_trajectory_time_cluster.svg", width = 10, height = 10)
plot_cell_trajectory(epi_cds, color_by = "Pseudotime", show_branch_points = F) + 
        scale_color_colormap(colormap = colormaps$density, reverse = T) +
        facet_wrap(~cluster, nrow = 2)
dev.off()

svglite(file = "./figures/epi_trajectory_time_group.svg", width = 10, height = 6)
plot_cell_trajectory(epi_cds, color_by = "Pseudotime", show_branch_points = F) + 
        scale_color_colormap(colormap = colormaps$density, reverse = T) +
        facet_wrap(~group, nrow = 1)
dev.off()

