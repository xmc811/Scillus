
# These functions are used for data processing and analytics

load_10X_mito_cc <- function(dir = NULL, raw_data = NULL, gcol = 2, proj_name, min_cells = 5, org = "human", cc = T) {
        
        if (!is.null(dir))  raw_data <- Read10X(data.dir = dir, gene.column = gcol)
        
        data <- CreateSeuratObject(counts = raw_data, project = proj_name, min.cells = min_cells)
        
        mito_pattern <- ifelse(org == 'human', "^MT-", "^mt-")
        
        if (org == "mouse") {
                
                s.genes <- vlookup(cc.genes$s.genes, mm_hs, 2, 1)
                g2m.genes <- vlookup(cc.genes$g2m.genes, mm_hs, 2, 1)
                
                s.genes <- s.genes[!is.na(s.genes)]
                g2m.genes <- g2m.genes[!is.na(g2m.genes)]
                
        } else {
                
                s.genes <- cc.genes$s.genes
                g2m.genes <- cc.genes$g2m.genes
        }
        
        data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = mito_pattern)
        
        if (cc) data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
        
        data$group <- proj_name
        
        return(data)
}



filter_norm_10X <- function(dataset, nfeature = 500, mito = 10, nfeatures = 2000) {
        
        expr1 <- FetchData(dataset, vars = "nFeature_RNA")
        expr2 <- FetchData(dataset, vars = "percent.mt")
        
        dataset <- dataset[, which(x = expr1 > nfeature & expr2 < mito)]
        dataset %<>% 
                NormalizeData(verbose = FALSE) %>%
                FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures)
        
}

filter_sctrans_10X <- function(dataset, nfeature = 500, mito = 10, 
                               vars = c("percent.mt","nCount_RNA","S.Score","G2M.Score")) {
        
        expr1 <- FetchData(dataset, vars = "nFeature_RNA")
        expr2 <- FetchData(dataset, vars = "percent.mt")
        
        dataset <- dataset[, which(x = expr1 > nfeature & expr2 < mito)]
        dataset %<>% 
                SCTransform(vars.to.regress = vars, verbose = T)
        
}

find_doublets <- function(dataset, dims = 1:20, ratio = 0.05, resolution = 0.4, txt) {
        
        dataset %<>%
                ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"), verbose = T) %>%
                RunPCA(features = VariableFeatures(dataset)) %>%
                RunUMAP(dims = dims) %>%
                FindNeighbors(dims = dims) %>%
                FindClusters(resolution = resolution)
        
        ## pK Identification
        sweep.res <- paramSweep_v3(dataset, PCs = dims)
        sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        
        ## Homotypic Doublet Proportion Estimate
        homotypic.prop <- modelHomotypic(dataset$seurat_clusters)
        nExp_poi <- round(ratio*length(Idents(dataset)))
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        
        dataset <- doubletFinder_v3(dataset, PCs = dims, pN = 0.25, pK = 0.1, nExp = nExp_poi.adj, reuse.pANN = F)
        
        barcodes <- names(dataset@active.ident)[dataset[[paste("DF.classifications_0.25_0.1_", as.character(nExp_poi.adj), sep = "")]] == "Doublet"]
        
        paste(barcodes, txt, sep = "")
        
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


rename_cluster <- function(dataset, labels) {
        
        if (length(labels) != length(levels(Idents(dataset)))) {
                stop("Length of new names must be the same with old names.")
        } else {
                current.name <- levels(Idents(dataset))
                new.name <- labels
                
                Idents(dataset) <- plyr::mapvalues(x = Idents(dataset), from = current.name, to = new.name)
                return(dataset)
        }
}


find_diff_genes <- function(dataset, clusters, groups, logfc = 0.25) {
        
        dataset$celltype.group <- paste(Idents(object = dataset), dataset$group, sep = "_")
        Idents(object = dataset) <- "celltype.group"
        
        
        de <- list()
        
        for (i in seq(length(clusters))) {
                
                d <- FindMarkers(dataset, 
                                  ident.1 = paste(clusters[i], groups[2], sep = "_"),
                                  ident.2 = paste(clusters[i], groups[1], sep = "_"),
                                  logfc.threshold = logfc,
                                  assay = "RNA")
                d %<>%
                        add_column(feature = rownames(d), .before = 1) %>%
                        add_column(cluster = clusters[i], .after = 1)
                
                de[[i]] <- as_tibble(d)
        }
        
        de <- do.call("rbind", de)
        
        de

}

test_GSEA <- function(diff, clusters, pathway) {
        
        gsea_res <- list()
        
        for (i in seq(length(clusters))) {
                
                data <- diff %>%
                        filter(cluster == clusters[i]) %>%
                        right_join(mm_hs, by = c("feature" = "mouse")) %>%
                        replace_na(list(avg_logFC = 0)) %>%
                        distinct(human, .keep_all = T) %>%
                        arrange(desc(avg_logFC))
                
                l <- data$avg_logFC
                names(l) <- data$human
                
                res <- fgsea(pathways = pathway, l, minSize = 15, maxSize = 500, nperm = 100000)
                
                res %<>%
                        add_column(cluster = clusters[i], .before = 1)
                
                gsea_res[[i]] <- res
        }
        
        gsea_res <- do.call("rbind", gsea_res)
        
        return(as_tibble(gsea_res))
        
}




add_program_score <- function(dataset, features, org = "human", nbin = 20, ctrl = 10, name){
  
        if(org == "mouse"){
                prog_genes <- vlookup(features, mm_hs, 2, 1)
                prog_genes <- list(prog_genes[!is.na(prog_genes)])
        } else {
                ex_genes <- list(features)
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


seurat_to_monocle <- function(dataset, subset = F, clusters = NULL) {
        
        dataset$cluster <- Idents(dataset)
        
        data <- GetAssayData(dataset, assay = 'RNA', slot = 'counts')
        pd <- new('AnnotatedDataFrame', data = dataset@meta.data)
        fd <- new('AnnotatedDataFrame', 
                      data = data.frame(gene_short_name = row.names(data), row.names = rownames(data)))
        
        cds <- newCellDataSet(data,
                                  phenoData = pd,
                                  featureData = fd,
                                  lowerDetectionLimit = 0.5,
                                  expressionFamily = negbinomial.size())
        cds %<>% 
                estimateSizeFactors() %>%
                estimateDispersions() %>%
                detectGenes(min_expr = 0.1)
        
        if (subset) {
                cds <- cds[,row.names(subset(pData(cds), cluster %in% clusters))]
        }
        
        return(cds)

}

analyze_monocle <- function(cds, rev = F) {
        
        expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
        
        diff_test_res <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~group")
        
        ordering_genes <- row.names(subset(diff_test_res, qval < 10E-70))
        
        cds %<>% 
                setOrderingFilter(ordering_genes) %>%
                reduceDimension(max_components = 2, method = 'DDRTree') %>%
                orderCells(reverse = rev)
        
        return(cds)

}


# Test

