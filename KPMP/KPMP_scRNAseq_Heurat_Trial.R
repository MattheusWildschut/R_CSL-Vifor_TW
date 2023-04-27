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

# Explore data & load
# hfile = Connect("Input\\521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat")
# hfile$index()
# data = LoadH5Seurat("Input\\521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat")
data <- subset(data, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
# data.scaled = data@assays$RNA@scale.data
data.matrix = data@assays$RNA@data
data.plot = data.frame("Expression" = data.matrix["P2RY14",],
           "Celltype" = data@meta.data$subclass.l2,
           "Cellclass" = data@meta.data$class,
           "Disease" = data@meta.data$condition.l1,
           "Patient" = data@meta.data$patient) %>%
  group_by(Disease) %>% dplyr::mutate(Disease = paste0(Disease, " (", n_distinct(Patient), ")")) %>%
  group_by(Celltype) %>% dplyr::mutate(Celltype = paste0(Celltype, " (", length(Expression), ")")) 
data.sum = data.plot %>%
  group_by(Disease, Celltype, Cellclass) %>% dplyr::summarize(Expression = mean(Expression), Percentage = sum(Expression>0)/n())

data.sum = data.plot %>% group_by(Disease, Patient, Celltype) %>% 
  dplyr::summarize(Percentage = sum(Expression>0)/n(), Expression = mean(Expression)) %>%
  tibble %>% tidyr::complete(Patient, Celltype, fill = list("Expression" = 0, "Percentage" = 0)) %>%
  filter(Celltype == "cDC (12)")
  
data.num = data.plot %>%
  group_by(Patient) %>% 
  dplyr::summarize(AllCells = length(Expression),
                SelCells = length(Expression[Celltype == "cDC (12)"]))
ggplot(data.num, aes(x = Patient)) +
  geom_col(aes(y = AllCells), fill = "black", width = 0.5, position = position_nudge(x = -0.25)) +
  geom_col(aes(y = SelCells*max(AllCells)/max(SelCells)), fill = "grey80", width = 0.5, position = position_nudge(x = 0.25)) +
  scale_y_continuous(name = "All cells", sec.axis = sec_axis(~./max(data.num$AllCells)*max(data.num$SelCells), name = "cDC (12) cells"), expand = c(0,0)) +
  theme_bw() + theme(axis.title.y.right = element_text(color = "grey80"), axis.text.y.right = element_text(color = "grey80"), axis.ticks.y.right = element_line(color = "grey80"))

unique(data@meta.data$subclass.full)
unique(data@meta.data$subclass.l2)


install.packages("xlsx")
library(xlsx)

wb = createWorkbook()
addDataFrame(as.data.frame(data.num), sheet = createSheet(wb, "Sheet 1"), row.names=FALSE)
addDataFrame(as.data.frame(data.sum), sheet = createSheet(wb, "Sheet 2"), row.names=FALSE)
saveWorkbook(wb, "My_File2.xlsx")
  
df <- tibble(
  group = c(1:2, 1, 2),
  item_id = c(1:2, 2, 3),
  item_name = c("a", "a", "b", "b"),
  value1 = c(1, NA, 3, 4),
  value2 = 4:7
)
df %>% complete(group, item_id, item_name)

  
  
  
  
ggplot(data.plot, aes(x = Disease, y = Celltype)) +
  geom_point(aes(size = Percentage, col = Expression)) +
  scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(data.plot$Percentage))) +
  scale_color_viridis() +
  facet_grid(Cellclass ~ ., scales = "free_y", space = "free") +
  theme_bw()

data.plot = data.frame("Expression" = data.scaled["CD3E",], 
                       "Celltype" = data@meta.data$subclass.l2,
                       "Cellclass" = data@meta.data$class,
                       "Disease" = data@meta.data$condition.l1,
                       "Patient" = data@meta.data$patient) %>%
  group_by(Disease) %>% dplyr::mutate(Disease = paste0(Disease, " (", n_distinct(Patient), ")")) %>%
  group_by(Celltype) %>% dplyr::mutate(Celltype = paste0(Celltype, " (", length(Expression), ")"))
b = filter(data.plot, Celltype == "T (4590)" & Disease == "AKI (12)")
sum(b$Expression>0)/nrow(b)
data.sum = data.plot %>%
  group_by(Disease, Celltype, Cellclass, Patient) %>% dplyr::summarize(Percentage = sum(Expression>0)/n(), Expression = mean(Expression))

ggplot(data.plot, aes(x = Expression, y = Celltype)) +
  geom_violin(aes(fill = Disease)) +
  scale_x_continuous(expand = c(0,0)) +
  facet_grid(Cellclass ~ ., scales = "free_y", space = "free") +
  theme_bw()

unique(data@meta.data$condition.l2)

DimPlot(data, reduction = "ref.umap")
DotPlot(data, features = "PADI4")

a = data@assays$RNA@data["CD3E",data@meta.data$subclass.l2 == "T"]
sum(a > 0)/length(a)

max(data@meta.data$percent.mt)

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
