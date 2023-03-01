suppressMessages({
  library(iasvaExamples)
  library(SingleCellExperiment)
  library(splatter)
  library(scuttle)
  source("jive_speedup.R")
})

###

# Parameters that remain the same for all cell types
nGenes <- 5000
method = "single"
dropout.type = "experiment"
dropout.shape = -1
verbose = FALSE

###

tic("Data Simulation")
data_SCE <- splatSimulate(
  # Parameters to tweak in each simulation
  batchCells = c(500, 500),
  batch.facLoc = 0.2,
  batch.facScale = 0.1,
  dropout.mid = 0.05,
  # Parameters to remain the same
  nGenes = nGenes,
  method = "groups",
  seed = 1,
  dropout.type = dropout.type,
  dropout.shape = dropout.shape,
  group.prob = c(0.6, 0.3, 0.1),
  de.prob = c(0.01, 0.1, 0.5),
  de.downProb = c(0.01, 0.4, 0.9),
  de.facLoc = c(0.6, 1, 0.2),
  de.facScale = c(0.1, 0.4, 0.8),
  verbose = verbose
)
toc()

###
data_SCE <- logNormCounts(data_SCE)

batches <- data_SCE$Batch
# Frequency of batches in simulation
table(batches)

groups <- data_SCE$Group
# Frequency of cell types in simulation
table(groups)

#######################
# Choose count matrix #
#######################

# counts <- counts(data_SCE)
counts <- logcounts(data_SCE)

dim(counts)

### Batches
n_batches <- length(unique(data_SCE$Batch))
unq_batches <- unique(data_SCE$Batch)

jive_data <- NULL
all_groups <- NULL

for (i in 1:n_batches) {
  jive_data[[paste0("Batch", i)]] <- t(counts[, batches == unq_batches[i]])
  all_groups[[paste0("Batch", i)]] <- groups[batches == unq_batches[i]]
}

###

tic("JIVE v2 runtime")
JIVE_results <- jive_v2(jive_data, rankJ = 10, rankA = rep(15, length(jive_data)), method = "given", maxiter = 5000)
# JIVE v2 runtime: 26587.25 sec elapsed
toc()

summary(JIVE_results)
saveRDS(JIVE_results, file = "data/JIVE_v2_dataset1_j10_a15.rds")

###

tic("JIVE v2 runtime")
JIVE_results <- jive_v2(jive_data, rankJ = 5, rankA = rep(20, length(jive_data)), method = "given", maxiter = 5000)
# JIVE v2 runtime: 2918.72 sec elapsed
toc()

summary(JIVE_results)
saveRDS(JIVE_results, file = "data/JIVE_v2_dataset1_j5_a20.rds")
