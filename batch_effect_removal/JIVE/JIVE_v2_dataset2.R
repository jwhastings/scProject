suppressMessages({
  library(tictoc)
  library(scRNAseq)
  library(scater)
  library(Seurat)
  source("jive_speedup.R")
})

###

BacherTCellData <- BacherTCellData()
table(BacherTCellData$batch)
table(BacherTCellData$batch, BacherTCellData$new_cluster_names)

BacherTCellData <- BacherTCellData[, !(BacherTCellData$new_cluster_names %in% c("Cycling", "Type-1 IFN signature"))]
table(BacherTCellData$batch)
table(BacherTCellData$batch, BacherTCellData$new_cluster_names)

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
BacherTCellData <- BacherTCellData[, colData(BacherTCellData)$batch %in% c(14, 15)]

# Final SCE object
data_SCE <- BacherTCellData

###
batches <- ifelse(data_SCE$batch < 10, paste0("Batch0", data_SCE$batch), paste0("Batch", data_SCE$batch))
# Frequency of batches in simulation
table(batches)

clusters <- data_SCE$new_cluster_names
# Frequency of cell types in simulation
table(batches, clusters)

counts <- counts(data_SCE)

dim(counts)

### Batches
n_batches <- length(unique(batches))
unq_batches <- sort(unique(batches))

jive_data <- NULL
all_clusters <- NULL

for (i in 1:n_batches) {
  jive_data[[unq_batches[i]]] <- as.matrix(t(counts[, batches == unq_batches[i]]))
  all_clusters[[unq_batches[i]]] <- clusters[batches == unq_batches[i]]
}

CORES <- parallel::detectCores()

###

tic("JIVE v2 runtime (permutation)")
JIVE_results_perm <- jive_v2(jive_data, method = "perm", maxiter = 100000, CORES = CORES)
# JIVE v2 runtime (permutaiton): 
toc()

summary(JIVE_results_perm)

saveRDS(JIVE_results_perm, file = "data/JIVE_v2_dataset2_alt_perm.rds")


###

tic("JIVE v2 runtime (given)")
JIVE_results <- jive_v2(jive_data, rankJ = 5, rankA = rep(9, length(jive_data)), method = "given", maxiter = 10000, CORES = CORES)
# JIVE v2 runtime (given): 789.33 sec elapsed
toc()

summary(JIVE_results)

saveRDS(JIVE_results, file = "data/JIVE_v2_dataset2_b14_b15_j10_a20.rds")
