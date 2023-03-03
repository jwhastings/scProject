suppressMessages({
  library(iasvaExamples)
  library(SingleCellExperiment)
  source("jive_speedup.R")
})

###

data("Lawlor_Islet_scRNAseq_Read_Counts")
data("Lawlor_Islet_scRNAseq_Annotations")

anns <- Lawlor_Islet_scRNAseq_Annotations
anns_subset <- anns[which(anns$Cell_Type != "none"), ]

counts <- Lawlor_Islet_scRNAseq_Read_Counts
counts_subset <- counts[, which(anns$Cell_Type != "none")]

# row_variance <- rowVars(counts_subset)
# top_5k_genevar <- tail(sort(rowVars(counts_subset), index.return = T)$ix, 15000)
# counts_subset <- counts_subset[top_5k_genevar, ]

data_SCE <- SingleCellExperiment(list(counts = counts_subset))
colData(data_SCE)$Batch <- anns_subset$Batch
colData(data_SCE)$Cell_Type <- anns_subset$Cell_Type

###

batches <- data_SCE$Batch
# Frequency of batches in simulation
table(batches)

groups <- data_SCE$Cell_Type
# Frequency of cell types in simulation
table(groups)

rawcounts <- counts(data_SCE)
dim(rawcounts)

# Batch 1
b1 <- rawcounts[, batches == "B1"]
grp_b1 <- groups[batches == "B1"]
# Dimension of batch 1
dim(b1)
# Frequency of cell types in batch 1
table(grp_b1)

# Batch 2
b2 <- rawcounts[, batches == "B2"]
grp_b2 <- groups[batches == "B2"]
# Dimension of batch 2
dim(b2)
# Frequency of cell types in batch 2
table(grp_b2)

# Batch 3
b3 <- rawcounts[, batches == "B3"]
grp_b3 <- groups[batches == "B3"]
# Dimension of batch 3
dim(b3)
# Frequency of cell types in batch 3
table(grp_b3)

###

jive_data <- NULL
jive_data$Batch1 <- t(b1)
jive_data$Batch2 <- t(b2)
jive_data$Batch3 <- t(b3)

###

CORES <- parallel::detectCores()

tic("JIVE v2 runtime")
JIVE_results <- jive_v2(jive_data, rankJ = 6, rankA = rep(8, length(jive_data)), method = "given", maxiter = 5000, CORES = CORES)
# JIVE v2 runtime: 26587.25 sec elapsed
toc()

saveRDS(JIVE_results, file = "data/JIVE_v2_dataset2.rds")
