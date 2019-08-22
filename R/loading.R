# Loading data
# Mitochondrial and cell cycle scoring

stages <- c("Early","Advanced")
samples <- c("GFP","CD45")

sample_names <- c(outer(samples, stages, paste, sep = "_"))
sample_paths <- paste("./data/HPW", c("1A","1B","2A","2B"), "/", sep = "")

scRNA <- pmap(.l = list(sample_paths, rep(stages, each = 2)), .f = load_10X_mito_cc, raw_data = NULL, org = "mouse", gcol = 2, min_cells = 5)
names(scRNA) <- sample_names

plot_qc(scRNA, "nFeature_RNA") + scale_y_log10()
plot_qc(scRNA, "nCount_RNA") + scale_y_log10()
plot_qc(scRNA, "percent.mt")

# Filtering and Normalization

scRNA <- purrr::map(scRNA, filter_norm_10X)

# Integration

gfp_combined <- IntegrateData(anchorset = FindIntegrationAnchors(object.list = scRNA[c(1,3)], dims = 1:30), dims = 1:30)
cd45_combined <- IntegrateData(anchorset = FindIntegrationAnchors(object.list = scRNA[c(2,4)], dims = 1:30), dims = 1:30)

# Doublet analysis

doublets <- map2(.x = scRNA, .y = rep(c("_1", "_2"), each = 2), find_doublets, dims = 1:20, ratio = 0.05, resolution = 0.4)


# Wahl data processing

wahl <- read_tsv("./data/GSE111113.tsv")
wahl <- wahl[!duplicated(wahl$features),]
wahl <- as.data.frame(wahl)
wahl <- wahl[complete.cases(wahl),]
rownames(wahl) <- wahl$features
wahl$features <- NULL

wahl_stages <- unique(str_extract(colnames(wahl), "^[A-Za-z0-9]+"))
wahl_stages <- wahl_stages[1:4]

wahl_list <- list()

for (i in seq(length(wahl_stages))) {
        wahl_list[[i]] <- dplyr::select(wahl, matches(wahl_stages[i]))
}

scRNA_wahl <- pmap(.l = list(wahl_list, wahl_stages), .f = load_10X_mito_cc, dir = NULL, org = "mouse", gcol = 2, min_cells = 5, cc = F)
names(scRNA_wahl) <- wahl_stages

rm(wahl)
rm(wahl_list)

plot_qc(scRNA_wahl, "nFeature_RNA") + scale_y_log10()
plot_qc(scRNA_wahl, "nCount_RNA") + scale_y_log10()
plot_qc(scRNA_wahl, "percent.mt")

scRNA_wahl <- purrr::map(scRNA_wahl, filter_norm_10X)

# Test





