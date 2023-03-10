### .......

library(harmony)

# Final SCE object
data_SCE <- BacherTCellData
colData(data_SCE)$Batch <- factor(colData(data_SCE)$batch)
colData(data_SCE)$Cluster <- factor(colData(data_SCE)$new_cluster_names)
data_SCE <- logNormCounts(data_SCE)
data_SCE <- runPCA(data_SCE, ncomponents = 30, ntop = 2000)
data_SCE <- RunHarmony(data_SCE, "Batch")
data_SCE <- runTSNE(data_SCE, dimred = "HARMONY")
data_SCE <- runUMAP(data_SCE, dimred = "HARMONY")

integrated_p1 <- plot_grid(
  plotReducedDim(data_SCE, "HARMONY", color_by = "Batch") +
    scale_color_fivethirtyeight(name = "Batch"),
  plotReducedDim(data_SCE, "HARMONY", color_by = "Cluster") + scale_color_gdocs(name = "Cluster"),
  nrow=2, align = "v"
)

integrated_p2 <- plot_grid(
  plotTSNE(data_SCE, color_by = "Batch") +
    scale_color_fivethirtyeight(name = "Batch"),
  plotTSNE(data_SCE, color_by = "Cluster") + scale_color_gdocs(name = "Cluster"),
  nrow=2, align = "v"
)

integrated_p3 <- plot_grid(
  plotUMAP(data_SCE, color_by = "Batch") +
    scale_color_fivethirtyeight(name = "Batch"),
  plotUMAP(data_SCE, color_by = "Cluster") + scale_color_gdocs(name = "Cluster"),
  nrow=2, align = "v"
)
