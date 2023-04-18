setwd("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/AKI_Immune_scRNAseq")
source("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/SourceFile_TW.R")

library(Seurat)
library(SeuratDisk)
URINE = readRDS("Input/SO_all_urine_cells.rds")
RENAL = readRDS("Input/SO_kidney_urine_cells.rds")

DimPlot(URINE, reduction = "umap", group.by = "renalcellgroup", label = TRUE, repel = TRUE, label.size = 3) +
  DimPlot(URINE, reduction = "umap", group.by = "renalcelltype", label = TRUE, repel = FALSE, label.size = 3) +
  FeaturePlot(URINE, features = "P2RY14")
ggsave("UMAPs_P2RY14.pdf", width = 22, height = 6, units = "in")
VlnPlot(URINE, features = "P2RY14", group.by = "renalcelltype")
DotPlot(URINE, features = "P2RY14", group.by = "renalcellgroup", dot.min = 0.0001) +
  DotPlot(URINE, features = "P2RY14", group.by = "renalcelltype", dot.min = 0.0001)
ggsave("DotPlots_P2RY14.pdf", width = 8, height = 7, units = "in")

(DotPlot(URINE, features = "P2RY14", group.by = "renalcellgroup", dot.min = 0.0001) +
  DotPlot(URINE, features = "P2RY14", group.by = "renalcelltype", dot.min = 0.0001)) /
  (DotPlot(URINE, features = c("P2RY14", "CLEC4C", "IL3RA", "TCF4", "IRF7", "IRF8", "CBFA2T3", "BCL11A"),
        group.by = "renalcellgroup", dot.min = 0.0001) +
  DotPlot(URINE, features = c("P2RY14", "CLEC4C", "IL3RA", "TCF4", "IRF7", "IRF8", "CBFA2T3", "BCL11A"),
          group.by = "renalcelltype", dot.min = 0.0001))
ggsave("DotPlots_P2RY14_DC-CellMarkers.pdf", width = 10, height = 14, units = "in")

DotPlot(subset(URINE, renalcellgroup %in% c("Bcells", "Tcells", "MYEL")), features = "P2RY14", group.by = "renalcellgroup", split.by = "AKI_type",
        dot.min = 0.0001, cols = viridis(3)) +
  scale_size_continuous(limits = c(0,4)) +
  scale_color_viridis(discrete = TRUE)
a = URINE@meta.data$p

FeatureScatter(URINE, "P2RY14", "CLEC4C", jitter = TRUE)
