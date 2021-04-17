
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
#' @param measures A string vector -
#' @param return_df A logical - if \code{TRUE}, return the dataframe; if \code{FALSE}, return a list of the dataframe and valid measure names. Default value is \code{TRUE}.
#' 
#' @return A tibble
#' @export
#' 

get_measure_data <- function(dataset, 
                             measures,
                             return_df = TRUE) {
        
        meta_vars <- colnames(dataset@meta.data)
        
        var_genes <- rownames(dataset@assays$integrated@scale.data)
        all_genes <- dataset@assays$RNA@data@Dimnames[[1]]
        
        df <- dataset@reductions$umap@cell.embeddings %>%
                as.data.frame() %>%
                cbind(dataset@meta.data) %>%
                rownames_to_column(var = "barcode")
        
        idx <- c() 
        
        for (i in seq_along(1:length(measures))) {
                
                if (measures[i] %in% meta_vars) {
                        next
                } else if (measures[i] %in% var_genes) {
                        df[[measures[i]]] <- as.numeric(dataset@assays$integrated@scale.data[measures[i],])
                } else if (measures[i] %in% all_genes) {
                        df[[measures[i]]] <- as.numeric(dataset@assays$RNA@data[measures[i],])
                } else {
                        warning(paste("Measure name not found:", measures[i]),  
                                immediate. = TRUE)
                        idx <- append(idx, i)
                }
        }
        
        if(!is.null(idx)) measures <- measures[-idx]
        
        if (return_df) {
                return(df)
        } else {
                return(list(df, measures))
        }
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
                arrange(desc(.data$avg_log2FC)) %>%
                group_by(.data$cluster) %>%
                filter(row_number() <= n) %>%
                arrange(.data$cluster)
        
        return(df$gene)
}


#' Test if a vector of strings are valid colors
#' 
#' @param x A string vector.
#' @return A named vector of logicals.
#' @importFrom grDevices col2rgb
#' 

are_colors <- function(x) {
        sapply(x, function(X) {
                tryCatch(is.matrix(col2rgb(X)), 
                         error = function(e) FALSE)
        })
}

#' Setup the colors for plotting
#' 
#' @param pal A string or a string vector.
#' @param n An integer.
#' @return A string vector of colors.
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom purrr map2
#' 

set_colors <- function(pal, n) {
        
        if (all(pal %in% rownames(brewer.pal.info))) {
                num <- c()
                for (i in seq(length(pal))) {
                        num[i] <- brewer.pal.info[pal[i],][[1]]
                }
                full_pal <- do.call(c, map2(.x = num, .y = pal, .f = brewer.pal))
        } else if (all(are_colors(pal))) {
                full_pal <- pal
        } else {
                stop('Incorrect palette setup. Please input valid RColorBrewer palette names or color names.')
        }
                
        if (n <= length(full_pal)) {
                return(full_pal[1:n])
        } else {
                warning("Number of colors required exceeds palette capacity. RdYlBu spectrum will be used instead.", 
                        immediate. = TRUE)
                return(colorRampPalette(brewer.pal(11, "RdYlBu"))(n))
        }
}



