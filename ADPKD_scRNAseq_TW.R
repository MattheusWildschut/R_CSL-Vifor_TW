setwd("C:\\Users\\mwildschut\\OneDrive - Vifor Pharma AG\\Documents\\R\\Projects\\R_CSL-Vifor_TW")

library(Seurat)
library(SeuratDisk)

"GSE185948"

h5_files <- list.files(path = "Input//ADPKD_snRNAseq", pattern = "*.h5")
h5_read <- lapply(h5_files, Read10X_h5)
h5_seurat <- lapply(h5_read, CreateSeuratObject)
h5_seurat  <- merge(h5_seurat[[1]], y = h5_seurat[2:length(h5_seurat)], 
                    add.cell.ids = 1:length(h5_files), project = "project")
SaveH5Seurat(h5_seurat, "Input//ADPKD_snRNAseq//ADPKD_snRNAseq_Combined.h5Seurat", overwrite = TRUE)
data = LoadH5Seurat("Input//ADPKD_snRNAseq//ADPKD_snRNAseq_Combined.h5Seurat")
data2 = data %>% subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)

a = fread("Input//ADPKD_snRNAseq//GSE185948_series_matrix.txt.gz")

dataset = "GSE185948"
gse = getGEO(dataset, GSEMatrix = TRUE)[[1]]
