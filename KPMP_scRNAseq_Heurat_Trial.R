## KPMP scRNAseq data downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183276
## Data loading
setwd("C:\\Users\\mwildschut\\OneDrive - Vifor Pharma AG\\Documents\\R\\Projects\\R_CSL-Vifor_TW")
source("SourceFile_TW.R")

# install.packages("Seurat")
# remotes::install_github("mojaveazure/seurat-disk")
# install.packages("Matrix")

library(Seurat)
library(SeuratDisk)

data = LoadH5Seurat("Input\\GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat",
                    assays = list(RNA = "data"))
# data <- subset(data, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
# all.genes <- rownames(data)
# gc()
# rm(data)
# data <- ScaleData(data, features = all.genes)
# data <- ScaleData(LoadH5Seurat("Input\\GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat",
#                                assays = list(RNA = "data")), features = all.genes)
# saveRDS(data, "Input\\GSE183276_KPMP_Healthy-Injury_Scaled.rds")
# data = readRDS("Input\\GSE183276_KPMP_Healthy-Injury_Scaled.rds")


data.scaled = data@assays$RNA@scale.data
data.scaled = data@assays$RNA@data
data.plot = data.frame("Gene" = data.scaled["SLC9B2",], 
                       "Celltype" = data@meta.data$subclass.l1,
                       "Cellclass" = data@meta.data$class,
                       "Disease" = data@meta.data$condition.l1) %>%
  group_by(Disease, Celltype, Cellclass) %>% dplyr::summarize(Expression = mean(Gene), Percentage = sum(Gene>0)/n()) %>%
  mutate(Product = Expression*Percentage)

# DotPlot(data, features = "SLC9B2", split.by = "condition.l2")

ggplot(data.plot, aes(x = Disease, y = Celltype)) +
  geom_point(aes(size = Percentage, col = Expression)) +
  scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(data.plot$Percentage))) +
  scale_color_viridis(option = "E") +
  facet_grid(Cellclass ~ ., scales = "free_y", space = "free") +
  theme_bw()
unique(data@meta.data$condition.l2)

DimPlot(data, reduction = "ref.umap")

# DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

raw.data = readRDS("Input//GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Counts_03282022.RDS")
data = CreateSeuratObject(counts = raw.data, project = "KPMP", min.cells = 3, min.features = 200)
data@active.ident = factor(data@active.ident, labels = rep("All samples", length(levels(data@active.ident))))
data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(data), 10)
plot1 <- VariableFeaturePlot(data)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data, features = VariableFeatures(object = data))
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca")
DimHeatmap(data, dims = 1:5, cells = 500, balanced = TRUE)

ElbowPlot(data, ndims = 50)
data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:15)

data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = 0.5)
head(Idents(data), 5)

data <- RunUMAP(data, dims = 1:15)
DimPlot(data, reduction = "umap", label = TRUE)

VlnPlot(data, features = c("SLC9B2"), slot = "counts", log = TRUE)
FeaturePlot(data, features = c("SLC9B2"))

data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(data, features = top10$gene) + NoLegend()
