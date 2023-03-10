suppressMessages({
  library(tictoc)
  library(scRNAseq)
  library(scater)
  library(Seurat)
  library(tidyverse)
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
BacherTCellData <- BacherTCellData[, colData(BacherTCellData)$batch %in% c(18, 19)]

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
unq_batches <- unique(batches)

jive_data <- NULL
all_clusters <- NULL

for (i in 1:n_batches) {
  jive_data[[paste0("Batch", i)]] <- as.matrix(t(counts[, batches == unq_batches[i]]))
  all_clusters[[paste0("Batch", i)]] <- clusters[batches == unq_batches[i]]
}

CORES <- parallel::detectCores() / 2

###

for (joint in 1:10) {
  tic("JIVE v2 runtime (given)")
  JIVE_results <- jive_v2(jive_data, rankJ = joint, rankA = c(15, 14), method = "given", maxiter = 10000, CORES = CORES)
  # JIVE v2 runtime (given): 789.33 sec elapsed
  toc()
  
  summary(JIVE_results)
  
  saveRDS(JIVE_results, file = paste0("data/JIVE_v2_dataset2_b18_b19_j", joint,"_a_15_14.rds"))
}

var_expl <- data.frame()
for (joint in 1:10) {
  temp <- readRDS(file = paste0("data/JIVE_v2_dataset2_b18_b19_j", joint,"_a_15_14.rds"))
  
  var <- as.data.frame(summary(temp)$Variance) %>% rownames_to_column(var = "Structure") %>%
    mutate(
      JointRank = joint,
      Indiv1Rank = 15,
      Indiv2Rank = 14
    )
  
  var_expl <- bind_rows(var_expl, var)
}

write_csv(var_expl, "output/dataset2_jive_variance_explained.csv")

var_expl %>%
  pivot_longer(cols = c(Batch1, Batch2), names_to = "Batch", values_to = "VarExpl") %>%
  ggplot(aes(x = JointRank, y = VarExpl)) +
  geom_point(aes(color = Structure)) +
  facet_wrap(Batch ~ .)
