library(scRNAseq)
library(scater)
library(tidyverse)

alldat <- data.frame(listDatasets())

all_calls <-
  # alldat$Call[1:2]
  alldat$Call[!(
    alldat$Call %in% c(
      "StoeckiusHashingData(mode='mouse')",
      "StoeckiusHashingData(mode='human')",
      "StoeckiusHashingData(type='mixed')"
    )
  )]

all_metadata_list <- lapply(all_calls, function(x) {
  cat(x, "\n")
  data.frame(
    data = x,
    vars = dimnames(colData(eval(parse(text = x))))[[2]]
  )
})

all_metadata <- bind_rows(all_metadata_list)

write_csv(all_metadata, file = "scRNA_coldata.csv")

batch <- all_metadata[!is.na(str_extract(all_metadata$vars, regex("batch", ignore_case = T))), ]
tech <- all_metadata[!is.na(str_extract(all_metadata$vars, regex("tech", ignore_case = T))), ]

###

# Load data with batch or technology

BacherTCellData <- BacherTCellData()
CampbellBrainData <- CampbellBrainData()
GiladiHSCData <- GiladiHSCData(mode='rna')
LedergorMyelomaData <- LedergorMyelomaData()
KotliarovPBMCData <- KotliarovPBMCData()
MessmerESCData <- MessmerESCData()
PaulHSCData <- PaulHSCData()
ZilionisLungData <- ZilionisLungData('mouse')
WuKidneyData <- WuKidneyData()

AztekinTailData <- AztekinTailData()
AztekinTailData
colData(AztekinTailData)$batch_factor <- factor(colData(AztekinTailData)$batch)
table(batch = colData(AztekinTailData)$batch_factor)
AztekinTailData <- logNormCounts(AztekinTailData)

AztekinTailData <- runPCA(AztekinTailData, ncomponents = 30, ntop = 2000)
plotPCA(AztekinTailData, color_by = "batch_factor")
plotPCA(AztekinTailData, color_by = "cluster")

AztekinTailData <- runTSNE(AztekinTailData)
plotTSNE(AztekinTailData, color_by = "batch_factor")
plotTSNE(AztekinTailData, color_by = "cluster")

AztekinTailData <- runUMAP(AztekinTailData)
plotUMAP(AztekinTailData, color_by = "batch_factor")
plotUMAP(AztekinTailData, color_by = "cluster")

###

BacherTCellData
table(batch = colData(BacherTCellData)$batch)

CampbellBrainData
table(batch = colData(CampbellBrainData)$batches)

GiladiHSCData
table(batch = colData(GiladiHSCData)$Amp_batch_ID)
table(batch = colData(GiladiHSCData)$Seq_batch_ID)

LedergorMyelomaData
table(batch = colData(LedergorMyelomaData)$Amp_batch_ID)

KotliarovPBMCData
table(batch = colData(KotliarovPBMCData)$batch)

MessmerESCData
table(batch = colData(MessmerESCData)$`experiment batch`)

PaulHSCData
table(batch = colData(PaulHSCData)$Seq_batch_ID)
table(batch = colData(PaulHSCData)$Amp_batch_ID)
table(batch = colData(PaulHSCData)$Batch_desc)

ZilionisLungData
table(batch = colData(ZilionisLungData)$`Library prep batch`)

WuKidneyData
table(batch = colData(WuKidneyData)$Technology)
