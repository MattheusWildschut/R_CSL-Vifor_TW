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

## KPMP-Immune gene list plots ------------------------------------------------------------------------------------------------

# gene.list = c("P2RY14", "CLEC4C", "IL3RA", "TCF4", "IRF7", "IRF8", "CBFA2T3", "BCL11A")
gene.list = c("CCL17", "IL3RA", "LY96", "CXCR3", "PLG")

data.scRNA2 = subset(data.scRNA, subclass.l1 == "Immune")
data.scRNA3 = subset(data.scRNA2, subclass.l2 == "cycMNP")
a = as.data.frame(data.scRNA3@assays$RNA@data) %>% select(contains("LD")) %>% 
  rownames_to_column("Gene")
a["LY96",]

meta.data = data.scRNA2@meta.data %>%
  mutate(Patient = str_remove(orig.ident, "a$|b$|c$"),
         Celltype = factor(.[, "subclass.l2"]),
         Disease = .[, "diseasetype"])

data.scRNA.pct = Percent_Expressing(data.scRNA2, features = gene.list,
                                    group_by = "subclass.l2", split_by = "diseasetype") %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = 2:ncol(.), values_to = "Percentage") %>%
  mutate(Disease = str_replace(str_split(name, "_", simplify = TRUE)[,2], "\\.", " "), Celltype = str_split(name, "_", simplify = TRUE)[,1],
         Percentage = Percentage/100) %>%
  select(-name) %>%
  mutate(Celltype = as.character(factor(Celltype, levels = sort(unique(Celltype)), labels = sort(levels(meta.data$Celltype)))))

data.scRNA.discellcounts = meta.data %>% dplyr::count(Celltype, Disease, .drop = FALSE, name = "DiseaseCellCounts") %>%
  complete(Celltype, Disease, fill = list(DiseaseCellCounts = 0)) %>%
  mutate(Disease = factor(Disease, levels = c("LivingDonor", "CKD", "AKI"))) %>%
  arrange(Disease) %>%
  group_by(Celltype) %>% 
  dplyr::mutate(Celltype2 = paste0(Celltype, " (", paste(DiseaseCellCounts, collapse = "|"), ")"))
data.scRNA.cellcounts = meta.data %>% dplyr::count(Celltype, .drop = FALSE, name = "CellCounts") %>%
  left_join(data.scRNA.discellcounts, by = "Celltype", multiple = "all")
data.scRNA.patcounts = meta.data %>% group_by(Disease) %>% dplyr::summarise("Patients" = n_distinct(Patient))

data.scRNA.avg = AverageExpression(subset(data.scRNA2, features = gene.list),
                                   group.by = c("subclass.l2", "diseasetype"))[[1]]

data.disease = data.scRNA.avg %>%
  t %>% as.data.frame %>%
  data.frame(., "Disease" = str_split_fixed(colnames(data.scRNA.avg), "_", 2)[,2], "Celltype" = str_split_fixed(colnames(data.scRNA.avg), "_", 2)[,1]) %>%
  pivot_longer(cols = 1:length(gene.list), names_to = "Gene", values_to = "Expression") %>%
  full_join(data.scRNA.pct, by = c("Gene", "Disease", "Celltype")) %>%
  full_join(data.scRNA.cellcounts, by = c("Disease", "Celltype")) %>%
  full_join(data.scRNA.patcounts, by = "Disease") %>%
  mutate(ClusterClass = capitalize(meta.data$ClusterClass[match(Celltype, meta.data$subclass.l2)]),
         Disease = paste0(Disease, " (", Patients, ")"),
         Celltype3 = paste0(Celltype, " (", CellCounts, ")")) %>%
  mutate(Expression = ifelse(DiseaseCellCounts < 10, NA, Expression))

data.counts = data.disease %>% 
  select(Disease, Celltype, CellCounts, DiseaseCellCounts) %>%
  distinct() %>%
  group_by(Disease) %>% 
  dplyr::mutate(TotalCells = sum(DiseaseCellCounts),
         Percentage = DiseaseCellCounts/TotalCells)
ggplot(data.counts, aes(y = Celltype)) +
  # scale_fill_viridis() +
  # new_scale_fill() +
  geom_tile(aes(x = Disease, fill = Percentage), col = "grey40") +
  scale_fill_distiller(direction = 1) +
  geom_text(aes(x = Disease, label = DiseaseCellCounts), 
            col = ifelse(data.counts$Percentage > quantile(data.counts$Percentage, 0.95), "white", "black")) +
  
  # scale_fill_viridis(option = "F", labels = scales::percent) +
  theme_bw() + scale_y_discrete(limits = rev) + coord_cartesian(expand = 0)


plot.list = map(as.list(gene.list), function(gene){
  data.gene = filter(data.disease, Gene == gene)
  ggplot(data.gene, aes(x = Disease, y = Celltype2)) +
    geom_point(aes(size = Percentage, col = Expression)) +
    # geom_text(aes(label = DiseaseCellCounts), size = 2, col = "grey50") +
    scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(c(0, data.gene$Percentage), na.rm = TRUE)), guide = guide_legend(direction = "vertical")) +
    scale_color_viridis(guide = guide_colorbar(direction = "vertical", order = 1), na.value="grey90") +
    scale_x_discrete(limits = rev) + scale_y_discrete(limits = rev) +
    facet_grid(ClusterClass ~ Gene, scales = "free", space = "free", switch = "y") + xlab(NULL) + ylab(NULL) +
    theme_bw() + theme(legend.position = "bottom", legend.box = "vertical",   #text = element_text(size = 18), 
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    {if(gene != gene.list[1]) labs(y = NULL)} +
    {if(gene != gene.list[1]) theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                    strip.background.y = element_blank(), strip.text.y.left = element_blank())}
})
wrap_plots(plot.list, nrow = 1)
ggsave("Output/KPMPquery-Immune-scRNA_CSLassets_TW.pdf", width = 9, height = 9)

## KPMP ClinicaL variables -----------------------------------------------------------------
raw.data = readxl::read_excel("Input/Biomarkers/OP-BC AU5812-Plasma Biomarker Data File-2022.xlsx",
                 sheet = 1)
raw.data2 = readxl::read_excel("Input/Biomarkers/OP-MSDQ120-Plasma Biomarker Data-2022.xlsx",
                              sheet = 1)
data = raw.data %>%
  `colnames<-`(str_remove_all(colnames(.), "\\ ")) %>%
  filter(`Sample Type` == "Serum") %>%
  left_join(raw.data2, by = "Participant ID") %>%
  # select(Participant.ID., Tissue.Type., Sample.Type., "Calcium..mg.dL..", "Phosphate..mg.dL..",
         # "Creatinine..mg.dL..") %>%
  mutate(`Calcium (mg/dL)` = as.numeric(`Calcium (mg/dL)`),
         `Phosphate (mg/dL)` = as.numeric(`Phosphate (mg/dL)`),
         Disease = factor(`Tissue Type`, levels = c("Healthy Reference", "AKI", "CKD")),
         Calcium_x_Phosphate = `Calcium (mg/dL)`*`Phosphate (mg/dL)`)
vars = list("Calcium (mg/dL)", "Phosphate (mg/dL)", "Calcium_x_Phosphate")

plot.list = map(vars, function(var){
  ggplot(data, aes(x = Disease, y = .data[[var]], col = Disease, fill = Disease)) +
    geom_boxplot(alpha = 0.25, outlier.alpha = 0) + geom_jitter(alpha = 0.5, width = 0.25) +
    stat_compare_means(comparisons = list(c(1,2), c(1,3), c(2,3))) +
    theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + xlab(NULL)
})
wrap_plots(plot.list, guides = "collect")

ggsave("Output/KPMP_Calcium-Phosphate_TW.pdf", width = 9, height = 3.5)

vars = list("Creatinine (mg/dL)", "Cystatin C (mg/L)",
            "C-reactive protein (high sensitivity) (mg/L )", "IL-6 Concentration (pg/mL)")
plot.list = map(vars, function(var){
  ggplot(data, aes(x = Calcium_x_Phosphate, y = .data[[var]], col = Disease)) +
    geom_point() +
    geom_smooth(formula = y~x, method = "lm") +
    stat_cor() +
    theme_bw()
})
wrap_plots(plot.list, guides = "collect")
data$`IL-6 Concentration (pg/mL)`
