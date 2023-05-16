setwd("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/KPMP")
source("../SourceFile_TW.R")

# install.packages("Seurat")
# remotes::install_github("mojaveazure/seurat-disk")
# install.packages("Matrix")

annot.scRNA = fread("Input/GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_Field_Descriptions.txt.gz") 
annot.snRNA = fread("Input/GSE183277_Kidney_Healthy-Injury_Cell_Atlas_snCv3_Metadata_Field_Descriptions.txt.gz")
pat.annot2 = read.csv("Input/72ec8f42-750f-4e1b-bda9-7300d4de1365_20221208_OpenAccessClinicalData.csv")

## Data from https://atlas.kpmp.org/repository/?facetTab=files&files_size=60&files_sort=%5B%7B%22field%22%3A%22
## file_name%22%2C%22order%22%3A%22asc%22%7D%5D&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in
## %22%2C%22content%22%3A%7B%22field%22%3A%22dois%22%2C%22value%22%3A%5B%2210.48698%2F92nk-e805%22%5D%7D%7D%5D%7D
data.scRNA = LoadH5Seurat("Input/521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat", assays = list(RNA = c("data", "counts")))
data.scRNA@meta.data$class = data.scRNA@meta.data$ClusterClass
data.scRNA@meta.data$Patient = str_remove(data.scRNA@meta.data$orig.ident, "a$|b$|c$")



## KPMP single cell correlations -------------------------------------------------------------------------------------------------
data.sel = subset(data.scRNA, subclass.l2 == "pDC")# & P2RY14 > 0)
gene = "CLEC4C"
scatter.plots = map(list("BCL11A", "CBFA2T3", "CLEC4C", "IRF7", "IRF8", "IL3RA", "TCF4", "CD83", "CD86"), function(gene){
  # data.sel
  FeatureScatter(data.sel, #cells = colnames(data.sel@assays$RNA@data)[data.sel@assays$RNA@data[gene,] > 0],
                 feature1 = "P2RY14", feature2 = gene, plot.cor = FALSE, jitter = TRUE)# +
    # geom_smooth(formula = y ~ x, method = "lm") +
    # stat_cor()
})
wrap_plots(scatter.plots, guides = "collect", nrow = 2)
ggsave("P2Y14/Output/P2RY14_pDC-Markers_CorrelationPlots_AllCells.pdf", width = 14, height = 6)



## KPMP gene list plots ------------------------------------------------------------------------------------------------

# gene.list = c("P2RY14", "CLEC4C", "IL3RA", "TCF4", "IRF7", "IRF8", "CBFA2T3", "BCL11A")
gene.list = c("CCL17", "IL3RA", "LY96", "CXCR3", "PLG")

meta.data = data.scRNA@meta.data %>%
  mutate(Patient = str_remove(orig.ident, "a$|b$|c$"),
         Celltype = factor(.[, "subclass.l1"]),
         Disease = .[, "diseasetype"])

data.scRNA.pct = Percent_Expressing(data.scRNA, features = gene.list,
                                    group_by = "subclass.l1", split_by = "diseasetype") %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = 2:ncol(.), values_to = "Percentage") %>%
  mutate(Disease = str_replace(str_split(name, "_", simplify = TRUE)[,2], "\\.", " "), Celltype = str_split(name, "_", simplify = TRUE)[,1],
         Percentage = Percentage/100) %>%
  select(-name) %>%
  mutate(Celltype = as.character(factor(Celltype, levels = sort(unique(Celltype)), labels = sort(levels(meta.data$Celltype)))))

data.scRNA.cellcounts = meta.data %>% dplyr::count(Celltype, .drop = FALSE, name = "CellCounts")
data.scRNA.patcounts = meta.data %>% group_by(Disease) %>% dplyr::summarise("Patients" = n_distinct(Patient))

data.scRNA.avg = AverageExpression(subset(data.scRNA, features = gene.list),
                                   group.by = c("subclass.l1", "diseasetype"))[[1]]

data.disease = data.scRNA.avg %>%
  t %>% as.data.frame %>%
  data.frame(., "Disease" = str_split_fixed(colnames(data.scRNA.avg), "_", 2)[,2], "Celltype" = str_split_fixed(colnames(data.scRNA.avg), "_", 2)[,1]) %>%
  pivot_longer(cols = 1:length(gene.list), names_to = "Gene", values_to = "Expression") %>%
  full_join(data.scRNA.pct, by = c("Gene", "Disease", "Celltype")) %>%
  full_join(data.scRNA.cellcounts, by = "Celltype") %>%
  full_join(data.scRNA.patcounts, by = "Disease") %>%
  mutate(ClusterClass = meta.data$ClusterClass[match(Celltype, meta.data$subclass.l1)],
         Disease = paste0(Disease, " (", Patients, ")"),
         Celltype = paste0(Celltype, " (", CellCounts, ")"))

plot.list = map(as.list(gene.list), function(gene){
  data.gene = filter(data.disease, Gene == gene)
  ggplot(data.gene, aes(x = Disease, y = Celltype)) +
    geom_point(aes(size = Percentage, col = Expression)) +
    scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(c(0, data.gene$Percentage), na.rm = TRUE)), guide = guide_legend(direction = "vertical")) +
    scale_color_viridis(guide = guide_colorbar(direction = "vertical", order = 1)) +
    scale_x_discrete(limits = rev) + scale_y_discrete(limits = rev) +
    facet_grid(ClusterClass ~ Gene, scales = "free", space = "free", switch = "y") + xlab(NULL) + ylab(NULL) +
    theme_bw() + theme(legend.position = "bottom", legend.box = "vertical",   #text = element_text(size = 18), 
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text.y.left = element_text(angle = 0)) +
    {if(gene != gene.list[1]) labs(y = NULL)} +
    {if(gene != gene.list[1]) theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                    strip.background.y = element_blank(), strip.text.y.left = element_blank())}
})
wrap_plots(plot.list, nrow = 1)
ggsave("Output/KPMPquery-scRNA_CSLassets_TW.pdf", width = 9, height = 9)
