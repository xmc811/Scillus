
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