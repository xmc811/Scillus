
# Helper functions

vlookup <- function(list, table, col_1, col_2) {
        new_list <- table[[col_2]][match(list, table[[col_1]])]
        return(new_list)
}




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



get_palette <- function(ncolor, palette = c("Paired", "Set2", "Set1")) {
        
        num <- c()
        for (i in seq(length(palette))) {
                num[i] <- brewer.pal.info[palette[i],][[1]]
        }
        
        ful_pal <- do.call(c, map2(.x = num, .y = palette, .f = brewer.pal))
        
        pal <- ful_pal[1:ncolor]
        
        return(pal)
        
}

#' Get a vector of colors.
#' 
#' @param v An integer index vector.
#' @param pal A string value describing the name of palette in \code{RColorBrewer}.
#' @return A vector of colors.
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @export
#' @examples
#' get_colors(c(1:2,4,6), "Paired")


get_colors <- function(v, pal = "Paired") {
        
        ncol <- brewer.pal.info[pal,][[1]]
        
        if (sum(!(v %in% 1:ncol)) > 0) {
                stop("Please input a valid numeric vector")
        }
        return(brewer.pal(n = ncol, name = pal)[v])
        
}
