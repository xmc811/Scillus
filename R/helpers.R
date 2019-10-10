
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
                filter(gene %in% int_features) %>%
                arrange(desc(avg_logFC)) %>%
                group_by(cluster) %>%
                filter(row_number() <= n) %>%
                arrange(cluster)
        
        return(df$gene)
}

#' Get a vector of colors for cluster coding
#' 
#' @param ncolor Number of colors to be generated
#' @param palette A string vector describing the name of palette in \code{RColorBrewer}.
#' @return A vector of colors.
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom purrr map2
#' @importFrom grDevices colorRampPalette
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
#' @param n Number of colors to be generated.
#' @return A vector of colors.
#' @export
#' 

get_spectrum <- function(n) {
        
        colorRampPalette(brewer.pal(12, "Set3"))(n)

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


