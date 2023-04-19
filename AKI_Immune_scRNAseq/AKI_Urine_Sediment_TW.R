setwd("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/AKI_Immune_scRNAseq")
source("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/SourceFile_TW.R")

library(Seurat)
library(SeuratDisk)
library(scCustomize)

URINE = readRDS("Input/SO_all_urine_cells.rds")
# RENAL = readRDS("Input/SO_kidney_urine_cells.rds")

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
  (DotPlot(URINE, features = c("P2RY14", "CLEC4C", "IL3RA", "TCF4", "IRF7", "IRF8", "CBFA2T3", "BCL11A", "ATP6V1C2"),
        group.by = "renalcellgroup", dot.min = 0.0001) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  DotPlot(URINE, features = c("P2RY14", "CLEC4C", "IL3RA", "TCF4", "IRF7", "IRF8", "CBFA2T3", "BCL11A", "ATP6V1C2"),
          group.by = "renalcelltype", dot.min = 0.0001) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
ggsave("DotPlots_P2RY14_DC-CellMarkers.pdf", width = 10, height = 14, units = "in")

DotPlot(subset(URINE, renalcellgroup %in% c("Bcells", "Tcells", "MYEL")), features = "P2RY14", group.by = "renalcellgroup", split.by = "AKI_type",
        dot.min = 0.0001) +
  scale_size_continuous(limits = c(0,4))# +
  # scale_color_viridis(discrete = TRUE)
a = URINE@meta.data$p

FeatureScatter(URINE, "P2RY14", "CLEC4C", jitter = TRUE)

## New trial with KPMP Shiny quantification script ----------------------------------
gene.list = c("P2RY14", "CLEC4C", "IL3RA", "TCF4", "IRF7", "IRF8", "CBFA2T3", "BCL11A")
pat.grouping = "AKI_type"
cell.grouping = "renalcellgroup"
scaling.list = FALSE
reorder.list = TRUE

data = URINE
data@meta.data = data@meta.data %>%
  mutate(Celltype = .[,cell.grouping],
         Disease = .[, pat.grouping],
         Patient = str_extract(patient, "P[:digit:]{3}|POOL[:digit:]"))

plot.list = map(list("MYEL", "Bcells", "Tcells"), function(cell){
  data.pct.pat = Percent_Expressing(data, features = gene.list, assay = "SCT",
                                    group_by = "Patient", split_by = cell.grouping) %>%
    `colnames<-`(str_replace(str_remove(colnames(.), "^X"), "\\.", "-")) %>%
    rownames_to_column("Gene") %>%
    pivot_longer(cols = 2:ncol(.), names_to = "Patient", values_to = "Percentage") %>%
    mutate(Percentage = Percentage/100,
           Celltype = str_replace(str_remove(Patient, "P[:digit:]{3}_|POOL[:digit:]_"), "CNT-CD", "CNT/CD"),
           Patient = str_extract(Patient, "P[:digit:]{3}|POOL[:digit:]"))
  
  data.counts.pat = data@meta.data %>%
    dplyr::count(Patient, Celltype, .drop = FALSE, name = "CellCounts") %>%
    left_join(dplyr::count(data@meta.data, Patient, name = "AllCells"), by = "Patient") %>%
    left_join(data@meta.data %>% group_by(Disease) %>% dplyr::mutate("Patients" = n_distinct(Patient)) %>%
                dplyr::select(Patient, Patients) %>% distinct, by = "Patient") %>%
    full_join(data.frame("Gene" = rep(gene.list, each = length(unique(data@meta.data$Patient))),
                         "Patient" = unique(data@meta.data$Patient)), by = "Patient", multiple = "all")
  
  data.avg.pat = AverageExpression(subset(data, features = gene.list), assay = "SCT",
                                   group.by = c("Patient", "Celltype"),
                                   slot = ifelse(scaling.list, "scale.data", "data"))[[1]]
  
  data.patient = data.avg.pat %>%
    t %>% as.data.frame %>% `colnames<-`(gene.list) %>%
    mutate(Patient = str_extract(rownames(.), "P[:digit:]{3}|POOL[:digit:]"),
           Celltype = str_remove(rownames(.), "P[:digit:]{3}_|POOL[:digit:]_")) %>%
    pivot_longer(cols = 1:(length(gene.list)), names_to = "Gene", values_to = "Expression") %>%
    full_join(data.pct.pat, by = c("Celltype", "Gene", "Patient")) %>%
    full_join(data.counts.pat, by = c("Celltype", "Gene", "Patient")) %>%
    mutate(Patient = paste0(Patient, " (", CellCounts, "/", AllCells, ")"),
      Disease = paste0(Disease, " (", Patients, ")"),
      Expression = ifelse(is.infinite(Expression), NA, Expression)) %>%
    filter(Celltype == cell)
  
  plot.patient = ggplot(data.patient, aes(x = if(reorder.list) reorder(Patient, -CellCounts) else Patient, y = Gene)) +
    geom_point(aes(size = Percentage, col = Expression)) +
    facet_grid(Celltype~Disease, scales = "free", space = "free", switch = "x") +
    scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(c(0, data.patient$Percentage), na.rm = TRUE)), guide = guide_legend(direction = "vertical")) +
    scale_color_viridis(guide = guide_colorbar(direction = "vertical")) + #, limits = quantile(data.patient$Expression, c(0.01,0.99), na.rm = TRUE)) +
    theme_bw() + theme(text = element_text(size = 18), strip.text.y.left = element_text(angle = 0), legend.position = "bottom",
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("Patient")
  
  data.pct.dis = Percent_Expressing(data, features = gene.list, assay = "SCT",
                                    group_by = pat.grouping, split_by = cell.grouping) %>%
    `colnames<-`(str_replace(str_remove(colnames(.), "^X"), "\\.", "-")) %>%
    rownames_to_column("Gene") %>%
    pivot_longer(cols = 2:ncol(.), names_to = "Disease", values_to = "Percentage") %>%
    mutate(Percentage = Percentage/100,
           Celltype = str_replace(str_remove(Disease, "CS_|pneumonia_|prerenal_"), "CNT-CD", "CNT/CD"),
           Disease = str_extract(Disease, "CS|pneumonia|prerenal"))
  
  data.counts.dis = data@meta.data %>%
    dplyr::count(Disease, Celltype, .drop = FALSE, name = "CellCounts") %>%
    left_join(dplyr::count(data@meta.data, Disease, name = "AllCells"), by = "Disease") %>%
    left_join(data@meta.data %>% group_by(Disease) %>% dplyr::mutate("Patients" = n_distinct(Patient)) %>%
                dplyr::select(Disease, Patients) %>% distinct, by = "Disease") %>%
    full_join(data.frame("Gene" = rep(gene.list, each = length(unique(data@meta.data$Disease))),
                         "Disease" = unique(data@meta.data$Disease)), by = "Disease", multiple = "all")
  
  data.avg.dis = AverageExpression(subset(data, features = gene.list), assay = "SCT",
                                   group.by = c("Disease", "Celltype"),
                                   slot = ifelse(scaling.list, "scale.data", "data"))[[1]]
  
  data.disease = data.avg.dis %>%
    t %>% as.data.frame %>% `colnames<-`(gene.list) %>%
    mutate(Disease = str_extract(rownames(.), "CS|pneumonia|prerenal"),
           Celltype = str_remove(rownames(.), "CS_|pneumonia_|prerenal_")) %>%
    pivot_longer(cols = 1:(length(gene.list)), names_to = "Gene", values_to = "Expression") %>%
    full_join(data.pct.dis, by = c("Celltype", "Gene", "Disease")) %>%
    full_join(data.counts.dis, by = c("Celltype", "Gene", "Disease")) %>%
    mutate(Disease = paste0(Disease, " (", CellCounts, "/", AllCells, ")"),
      Expression = ifelse(is.infinite(Expression), NA, Expression)) %>%
    filter(Celltype == cell) #"Bcells", "Tcells", 
  
  plot.disease = ggplot(data.disease, aes(x = if(reorder.list) reorder(Disease, -CellCounts) else Disease, y = Gene)) +
    geom_point(aes(size = Percentage, col = Expression)) +
    scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(c(0, data.disease$Percentage), na.rm = TRUE)), guide = guide_legend(direction = "vertical")) +
    scale_color_viridis(guide = guide_colorbar(direction = "vertical")) + #, limits = quantile(data.patient$Expression, c(0.01,0.99), na.rm = TRUE)) +
    theme_bw() + theme(text = element_text(size = 18), strip.text.y.left = element_text(angle = 0), legend.position = "bottom",
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("Disease")
  
  plot.disease + plot.patient + plot_layout(widths = c(1,7))
})
wrap_plots(plot.list, ncol = 1)

