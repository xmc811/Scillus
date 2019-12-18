# Thess functions are deprecated but may be useful
find_markers <- function(dataset, use.par = F, ncores = 4, min.diff.pct = -Inf) {
        
        DefaultAssay(object = dataset) <- "RNA"
        
        if (use.par == F) {
                
                markers <- list()
                
                for (i in 0:(length(levels(Idents(dataset)))-1)) {
                        
                        m <- FindConservedMarkers(object = dataset, ident.1 = i, grouping.var = "group", verbose = T, min.diff.pct = min.diff.pct)
                        m %<>%
                                add_column(feature = rownames(m), .before = 1) %>%
                                add_column(cluster = i, .after = 1)
                        
                        markers[[i+1]] <- m
                }
                
        } else {
                
                registerDoParallel(cores = ncores)
                
                markers <- foreach(i = 0:(length(levels(Idents(dataset)))-1)) %dopar% {
                        
                        m <- FindConservedMarkers(object = dataset, ident.1 = i, grouping.var = "group", verbose = T, min.diff.pct = min.diff.pct)
                        m %<>%
                                add_column(feature = rownames(m), .before = 1) %>%
                                add_column(cluster = i, .after = 1)
                        m
                }
        }
        markers <- do.call("rbind", markers)
        
        rownames(markers) <- NULL
        markers <- as_tibble(markers)
        
        return(markers)
}


analyze_merged <- function(dataset, group.levels,
                           verbose = T, npcs = 50,
                           reduction = "umap", dims = 1:20, nnei = 30, min.dist = 0.3, spread = 1, n.epochs = 500, 
                           k.param = 20,
                           resolution = 0.8) {
        
        dataset$group <- factor(dataset$group, levels = group.levels)
        dataset %<>% 
                RunPCA(npcs = npcs, verbose = verbose) %>%
                RunUMAP(reduction = "pca", dims = dims, n.neighbors = nnei, min.dist = min.dist, spread = spread, n.epochs = n.epochs) %>%
                FindNeighbors(reduction = reduction, dims = 1:2, k.param = k.param) %>%
                FindClusters(resolution = resolution, algorithm = 3)
        
}
