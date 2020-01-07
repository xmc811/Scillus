
# These functions are used for data processing and analytics

#' Load single cell RNA-seq data from directory or file
#' 
#' @param metadata A dataframe - must contain a column named "sample", and a column named "folder" or "file".
#' @param ... Additional arguments to be passed to the function \code{\link{CreateSeuratObject}}
#' 
#' @return A list of Seurat object
#' 
#' @importFrom Seurat Read10X CreateSeuratObject PercentageFeatureSet
#' @importFrom purrr map
#' @importFrom data.table fread
#' @export
#' 

load_scfile <- function(metadata, ...) {
        
        if (!"sample" %in% colnames(metadata)) {
                stop("No sample names in the metadata")
        }
        
        a <- "folder" %in% colnames(metadata) + "file" %in% colnames(metadata)
        
        if (a == 0) {
                stop("No file or folder names in the metadata")
        } else if (a > 1) {
                warning("Both file and folder names are provided in the metadata. Only folder names will be processed", 
                        immediate. = TRUE)
        }
        
        if (any(duplicated(metadata[['sample']]))) {
                warning("Sample names are not unique", 
                        immediate. = TRUE)
        }
        
        if ("folder" %in% colnames(metadata)) {
                mat_list <- purrr::map(.x = metadata[['folder']], 
                                       .f = Read10X)
        } else {
                mat_list <- purrr::map(.x = metadata[['file']], 
                                       .f = function(x) data.frame(data.table::fread(file = x), 
                                                                   row.names = T))
        }
        
        data <- list()
        
        meta_var <- colnames(metadata)
        meta_var <- meta_var[!meta_var %in% c("sample", "folder", "file")]
        
        for (i in seq_along(1:length(mat_list))) {
                
                data[[i]] <- CreateSeuratObject(counts = mat_list[[i]], 
                                                project = as.character(metadata[['sample']][i]), ...)
                data[[i]][['percent.mt']] <- PercentageFeatureSet(data[[i]], 
                                                                  pattern = "(^(hg19|mm10)-(MT|mt)-|^(MT|mt)-)")
                data[[i]][['sample']] <- metadata[['sample']][i]
                
                if (length(meta_var) > 0) {
                        
                        for (j in seq_along(1:length(meta_var))) {
                                data[[i]][[meta_var[j]]] <- metadata[[meta_var[j]]][i]
                        }
                }
        }
        return(data)
}

#' Filter a list of Seurat objects
#' 
#' @param data_list A list of Seurat objects
#' @param ... Additional arguments to be passed to the function \code{\link{subset.Seurat}}
#' 
#' @return A list of Seurat objects
#' 
#' @importFrom ggplot2 position_dodge
#' @export
#' 

filter_scdata <- function(data_list, ...) {
        
        sample <- purrr::map_chr(.x = data_list, 
                                 .f = function(x) x@project.name)
        pre <- purrr::map_int(.x = data_list, 
                              .f = function(x) nrow(x[["nCount_RNA"]]))
        data_list <- purrr::map(.x = data_list, 
                                .f = subset, ...)
        post <- purrr::map_int(.x = data_list, 
                               .f = function(x) nrow(x[["nCount_RNA"]]))
        
        df <- tibble::tibble(sample = sample,
                             pre = pre,
                             post = post)
        
        df %<>%
                tidyr::gather(variable, value, -sample) %>%
                mutate(variable = factor(.data$variable, levels = c("pre","post")))
        
        p <- ggplot(df) +
                geom_bar(aes(x = sample,
                             y = .data$value,
                             fill = .data$variable), 
                         stat = "identity", 
                         position = "dodge") +
                scale_fill_manual(values = c("#fbb4ae","#b3cde3"), 
                                  labels = c("Pre-QC", "Post-QC"), 
                                  name = NULL) +
                geom_text(aes(x = sample, 
                              y = .data$value, 
                              label = .data$value, 
                              group = .data$variable), 
                          position = position_dodge(width = 1), 
                          vjust = -0.5, 
                          size = 3.5) +
                labs(y = "Number of Cells") +
                theme(legend.title = element_text(size = 12),
                      legend.text = element_text(size = 12),
                      axis.title = element_text(size = 12),
                      axis.text = element_text(size = 12),
                      axis.title.x = element_blank())
        
        graphics::plot(p)
        return(data_list)
}

#' Refactor grouping variables in Seurat object
#' 
#' @param dataset A Seurat object.
#' @param metadata A dataframe - should be the same as the one used during loading. Default value is \code{NULL}.
#' 
#' @return A Seurat object.
#' @export
#' 

refactor_seurat <- function(dataset, 
                            metadata = NULL) {
        
        if (is.null(metadata)) {
                
                message("No metadata provided. All the metadata of 'character' type will be factored.")
                fct_variable <- colnames(dataset@meta.data)[sapply(dataset@meta.data, is.character)]
                
                for (i in seq_along(1:length(fct_variable))) {
                        
                        dataset[[fct_variable[i]]][[1]] <- factor(dataset[[fct_variable[i]]][[1]])
                }
        } else {
                fct_variable <- colnames(metadata)[sapply(metadata, is.factor)]
                
                for (i in seq_along(1:length(fct_variable))) {
                        
                        dataset[[fct_variable[i]]][[1]] <- factor(dataset[[fct_variable[i]]][[1]], 
                                                                  levels = levels(metadata[[fct_variable[[i]]]]))
                }
        }
        return(dataset)
}



#' Rename clusters in Seurat object
#' 
#' @param dataset A Seurat object.
#' @param labels An string vector - new names of each cluster. 
#' 
#' @return A Seurat object.
#' @importFrom Seurat Idents
#' @importFrom plyr mapvalues
#' @export
#' 

rename_cluster <- function(dataset, labels) {
        
        if (length(labels) != length(levels(Idents(dataset)))) {
                stop("Length of new names must be the same with old names.")
        } else {
                current.name <- levels(Idents(dataset))
                new.name <- labels
                
                Seurat::Idents(dataset) <- plyr::mapvalues(x = Idents(dataset), from = current.name, to = new.name)
                return(dataset)
        }
}

#' Find differential expressed genes between two groups in each cluster
#' 
#' @param dataset A Seurat object.
#' @param clusters A string vector - clusters to investigate.
#' @param comparison A string vector of length 3 -
#' @param logfc A double -
#' 
#' @return A Seurat object.
#' @importFrom Seurat Idents FindMarkers
#' @importFrom tibble add_column as_tibble
#' @importFrom magrittr %>% %<>%
#' @export
#' 

find_diff_genes <- function(dataset, 
                            clusters, 
                            comparison, 
                            logfc = 0.25) {
        
        Seurat::Idents(dataset) <- paste(as.character(Idents(dataset)), dataset[[comparison[1]]][[1]], sep = "_")

        de <- list()
        
        for (i in seq_along(1:length(clusters))) {
                
                d <- FindMarkers(dataset, 
                                 ident.1 = paste(clusters[i], comparison[3], sep = "_"),
                                 ident.2 = paste(clusters[i], comparison[2], sep = "_"),
                                 logfc.threshold = logfc,
                                 assay = "RNA")
                de[[i]] <- d %>%
                        rownames_to_column(var = "feature") %>%
                        add_column(cluster = clusters[i], .after = 1) %>%
                        as_tibble()
        }
        return(do.call("rbind", de))
}

#' GSEA analysis of differential gene expression in each cluster
#' 
#' @param diff A tibble -
#' @param clusters A string vector - clusters to investigate. Default value is \code{NULL}.
#' @param pathway A vector list - group of gene lists 
#' 
#' @return A Seurat object.
#' @importFrom Seurat Idents FindMarkers
#' @importFrom dplyr filter right_join distinct arrange desc
#' @importFrom tidyr replace_na
#' @importFrom tibble add_column as_tibble
#' @importFrom magrittr %>% %<>%
#' @export
#' 

test_GSEA <- function(diff, 
                      clusters = NULL, 
                      pathway) {
        
        if (!requireNamespace("fgsea", quietly = TRUE)) {
                stop(paste("Package \"fgsea\" needed for this function to work. Please install it."),
                     call. = FALSE)
        }
        
        if (is.null(clusters)) clusters <- unique(diff$cluster)
        
        gsea_res <- list()
        
        for (i in seq(length(clusters))) {
                
                data <- diff %>%
                        filter(.data$cluster == clusters[i]) %>%
                        arrange(desc(.data$avg_logFC))
                
                l <- data$avg_logFC
                names(l) <- data$feature
                
                res <- fgsea::fgsea(pathways = pathway, 
                                    stats = l, 
                                    minSize = 15, 
                                    maxSize = 500, 
                                    nperm = 10000)
                
                res %<>%
                        add_column(cluster = clusters[i],
                                   .before = 1)
                
                gsea_res[[i]] <- res
        }
        return(as_tibble(do.call("rbind", gsea_res)))
}

#' Add gene program scores
#' 
#' @param dataset A Seurat object.
#' @param features An string vector - a gene list of expression programs. 
#' @param org A string - name of organism. Currently only "human" or "mouse" are accepted. 
#' @param nbin An integer - number of bins of aggregate expression levels for all analyzed features
#' @param ctrl An integer - number of control features selected from the same bin per analyzed feature
#' @param name A string - name of the expression program
#' 
#' @return A Seurat object.
#' @importFrom Seurat AddModuleScore
#' @export
#' 

add_program_score <- function(dataset, features, org = "human", nbin = 20, ctrl = 10, name){
        
        if(org == "mouse"){
                prog_genes <- vlookup(features, mm_hs, 2, 1)
                prog_genes <- list(prog_genes[!is.na(prog_genes)])
        } else {
                prog_genes <- list(features)
        }
        
        n_genes <- nrow(dataset@assays$integrated@scale.data)
        
        genes_per_bin <- round(n_genes/nbin)
        
        ctrl <- ifelse(ctrl > genes_per_bin, round(genes_per_bin/3), ctrl)
        
        dataset <- AddModuleScore(dataset,
                                  features = prog_genes,
                                  ctrl = ctrl,
                                  nbin = nbin,
                                  name = name)
        
        return(dataset)
}



