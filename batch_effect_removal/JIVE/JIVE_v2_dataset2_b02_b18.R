suppressMessages({
  library(tictoc)
  library(scRNAseq)
  library(scater)
  library(Seurat)
  source("jive_speedup.R")
})

###

BacherTCellData <- BacherTCellData()
bacher_metadata <- as.data.frame(colData(BacherTCellData))
bacher_seurat <- CreateSeuratObject(counts(BacherTCellData), meta.data = bacher_metadata)

#################
# Identify HVGs #
#################

# split the dataset into a list of seurat objects (one for each batch)
bacher_seurat_batch <- SplitObject(bacher_seurat, split.by = "batch")

# normalize and identify variable features for each dataset independently
bacher_seurat_batch <- lapply(X = bacher_seurat_batch, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = bacher_seurat_batch)

########################################################

# Select only HVGs
BacherTCellData <- BacherTCellData[features, ]

# Subset to three batches
BacherTCellData <- BacherTCellData[, colData(BacherTCellData)$batch %in% c(2, 18)]

# Final SCE object
data_SCE <- BacherTCellData

###
batches <- ifelse(data_SCE$batch < 10, paste0("Batch0", data_SCE$batch), paste0("Batch", data_SCE$batch))
# Frequency of batches in simulation
table(batches)

clusters <- data_SCE$new_cluster_names
# Frequency of cell types in simulation
table(clusters)

counts <- counts(data_SCE)

dim(counts)

### Batches
n_batches <- length(unique(batches))
unq_batches <- unique(batches)

jive_data <- NULL
all_clusters <- NULL

for (i in 1:n_batches) {
  jive_data[[paste0("Batch", i)]] <- t(counts[, batches == unq_batches[i]])
  all_clusters[[paste0("Batch", i)]] <- clusters[batches == unq_batches[i]]
}

CORES <- parallel::detectCores()

###

tic("JIVE v2 runtime (BIC)")
JIVE_results_bic <- jive_v2(jive_data, method = "bic", maxiter = 5000, CORES = CORES)
# JIVE v2 runtime: 6577.61 sec elapsed
toc()

saveRDS(JIVE_results, file = "data/JIVE_v2_dataset2_b02_b18_bic.rds")
