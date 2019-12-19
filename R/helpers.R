
# Helper functions

#' Look up corresponding values from a lookup table
#' 
#' @param list A string vector - keys
#' @param table A tibble - the lookup table
#' @param col_1 An integer - column index of the keys
#' @param col_2 An integer - column index of the values
#' 
#' @return A string vector.
#' @export
#' 

vlookup <- function(list, table, col_1, col_2) {
        new_list <- table[[col_2]][match(list, table[[col_1]])]
        return(new_list)
}


#' Conversion between mouse and human gene symbols
#' 
#' @param gene_list A string vector - a list of genes to convert
#' @param hs_to_mm A logical value - if the conversion is from human to mouse. Default value is \code{TRUE}.
#' 
#' @return A string vector.
#' @export
#' 

convert_symbol <- function(gene_list, hs_to_mm = TRUE) {
        
        genes <- vlookup(gene_list, mm_hs, 1 + hs_to_mm, 2 - hs_to_mm)
        genes[!is.na(genes)]
       
}


#' Test whether a gene waas in the dataset
#' 
#' @param dataset A Seurat object -
#' @param gene A string -
#' 
#' @return A double
#' @export
#' 

gene_in_data <- function(dataset, gene) {
        
        if (!(gene %in% dataset@assays$RNA@data@Dimnames[[1]])) {
                return(0)
        } else if (!(gene %in% dataset@assays$integrated@data@Dimnames[[1]])) {
                return(1)
        } else {
                return(2)
        }
}


#' Obtain the expression data for genes
#' 
#' @param dataset A Seurat object -
#' @param genes A string vector -
#' 
#' @return A tibble
#' @export
#' 


get_gene_data <- function(dataset, genes) {
        
        genes <- genes[map_dbl(.x = genes, .f = gene_in_data, dataset = dataset) > 0]
        
        lst <- list()
        
        for (i in seq_along(genes)) {
                
                if (gene_in_data(dataset, genes[i]) == 2) {
                        
                        lst[[i]] <- tibble(group = as.character(dataset$group),
                                           cluster = as.character(Idents(dataset)),
                                           x_cor = dataset@reductions$umap@cell.embeddings[,1],
                                           y_cor = dataset@reductions$umap@cell.embeddings[,2],
                                           value = dataset@assays$integrated@data[genes[i],],
                                           feature = genes[i])
                } else {
                        
                        lst[[i]] <- tibble(group = as.character(dataset$group),
                                           cluster = as.character(Idents(dataset)),
                                           x_cor = dataset@reductions$umap@cell.embeddings[,1],
                                           y_cor = dataset@reductions$umap@cell.embeddings[,2],
                                           value = dataset@assays$RNA@data[genes[i],],
                                           feature = genes[i])
                }
                
        }
        
        df <- do.call("rbind", lst)
        
        return(df)
}


#' Obtain the meta data for features
#' 
#' @param dataset A Seurat object -
#' @param measures A string vector -
#' 
#' @return A tibble
#' @export
#' 

get_meta_data <- function(dataset, measures) {
        
        measures <- measures[measures %in% colnames(dataset@meta.data)]
        
        lst <- list()
        
        for (i in seq_along(measures)) {
                        
                lst[[i]] <- tibble(group = as.character(dataset$group),
                                   cluster = as.character(Idents(dataset)),
                                   x_cor = dataset@reductions$umap@cell.embeddings[,1],
                                   y_cor = dataset@reductions$umap@cell.embeddings[,2],
                                   value = as.numeric(dataset@meta.data[[measures[i]]]),
                                   feature = measures[i])
                
        }
        
        df <- do.call("rbind", lst)
        
        return(df)
}



#' Get top genes of each cluster
#' 
#' @param dataset A Seurat object.
#' @param markers A tibble - from FindMarkers function. 
#' @param n An integer - number of top genes to choose. 
#' 
#' @return A string vector of genes.
#' @importFrom dplyr filter arrange group_by row_number
#' @importFrom magrittr %>% %<>%
#' @export
#' 

get_top_genes <- function(dataset, markers, n) {
        
        int_features <- rownames(dataset@assays$integrated@scale.data)
        
        df <- markers %>%
                filter(.data$gene %in% int_features) %>%
                arrange(desc(.data$avg_logFC)) %>%
                group_by(.data$cluster) %>%
                filter(row_number() <= n) %>%
                arrange(.data$cluster)
        
        return(df$gene)
}

#' Get a vector of colors for cluster coding
#' 
#' @param ncolor Number of colors to be generated
#' @param palette A string vector describing the name of palette in \code{RColorBrewer}.
#' @return A vector of colors.
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom purrr map2
#' @export
#' 
#' @examples
#' get_palette(12, c("Set2", "Set1"))
#' 

get_palette <- function(ncolor, palette = c("Paired", "Set2", "Set1")) {
        
        num <- c()
        for (i in seq(length(palette))) {
                num[i] <- brewer.pal.info[palette[i],][[1]]
        }
        
        ful_pal <- do.call(c, map2(.x = num, .y = palette, .f = brewer.pal))
        
        pal <- ful_pal[1:ncolor]
        
        return(pal)
        
}

#' Get a vector of colors for sample coding
#' 
#' @param n Number of colors to be generated
#' @param palette A string
#' @return A vector of colors
#' @importFrom grDevices colorRampPalette
#' @export
#' 

get_spectrum <- function(n, palette = "Set2") {
        
        if (n <= brewer.pal.info[palette,]['maxcolors'][[1]]) {
                return(suppressWarnings(brewer.pal(n, name = palette)))
        } else {
                warning("Number of colors exceeds palette capacity. RdYlBu spectrum will be used instead.", immediate. = TRUE)
                return(colorRampPalette(brewer.pal(11, "RdYlBu"))(n))
        }
}


#' Get a vector of colors
#' 
#' @param v An integer index vector.
#' @param pal A string value describing the name of palette in \code{RColorBrewer}.
#' @return A vector of colors.
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @export
#' 
#' @examples
#' get_colors(c(1:2,4,6), "Paired")
#' 

get_colors <- function(v, pal = "Paired") {
        
        ncol <- brewer.pal.info[pal,][[1]]
        
        if (sum(!(v %in% 1:ncol)) > 0) {
                stop("Please input a valid numeric vector")
        }
        return(brewer.pal(n = ncol, name = pal)[v])
        
}


