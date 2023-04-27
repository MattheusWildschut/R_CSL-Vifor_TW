## Data loading
setwd("C:\\Users\\mwildschut\\OneDrive - Vifor Pharma AG\\Documents\\R\\Projects\\R_CSL-Vifor_TW")
source("SourceFile_TW.R")

dataset = "GSE7869"
gse = getGEO(dataset, GSEMatrix = TRUE)[[1]]

uniprot = fread("Input\\Uniprot_HumanProteome_28-02-23.tsv")
annot = gse@featureData@data
raw.data = gse@assayData$exprs[annot$`Gene Symbol` != "" & !str_detect(annot$`Gene Symbol`, "-| ///|\\."),]
annot = annot[annot$`Gene Symbol` != "" & !str_detect(annot$`Gene Symbol`, "-| ///|\\."),]
annot.pat = gse@phenoData@data %>%
  mutate(Group = factor(case_when(str_detect(title, "normal|Normal") ~ "Normal (n=3)",
                                  str_detect(title, "minimally") ~ "MCD (n=5)",
                                  str_detect(title, "small") ~ "Cysts-Small (n=5)",
                                  str_detect(title, "medium") ~ "Cysts-Medium (n=5)",
                                  str_detect(title, "large") ~ "Cysts-Large (n=3)"),
                        levels = c("Normal (n=3)", "MCD (n=5)", "Cysts-Small (n=5)", "Cysts-Medium (n=5)", "Cysts-Large (n=3)")))

annot.genes = data.frame("Gene" = sort(unique(annot$`Gene Symbol`)), 
                         "Description" = annot$`Gene Title`[match(sort(unique(annot$`Gene Symbol`)), annot$`Gene Symbol`)])

max.data = as.data.frame(raw.data) %>%
  mutate(Gene = annot$`Gene Symbol`) %>%
  group_by(Gene) %>%
  dplyr::summarise(across(everything(), max)) %>%
  column_to_rownames("Gene") %>% as.matrix

data.cor = data.frame("Spearman" = cor(t(max.data), as.numeric(annot.pat$Group), method = "spearman"),
               "Pearson" = cor(t(max.data), as.numeric(annot.pat$Group), method = "pearson")) %>% 
  as.data.frame %>% rownames_to_column("Gene") %>% arrange(desc(Spearman)) %>%
  mutate("Rank_Spearman" = rank(-Spearman), "Rank_Pearson" = rank(-Pearson))

kegg2 = getKEGGPathwayNames(species="hsa")
kegg1 = getGeneKEGGLinks(species="hsa") %>%
  mutate(Gene = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, GeneID, column = "SYMBOL", keytype = "ENTREZID"),
         PathwayID = str_remove(PathwayID, "path:"),
         Description = str_remove(kegg2$Description[match(PathwayID, kegg2$PathwayID)], " - Homo sapiens \\(human\\)")) %>%
  select(Description, Gene)

GSEA.correlations = function(cor_method){
  list.kegg = setNames(data.cor[,cor_method], data.cor$Gene)[order(data.cor[,cor_method], decreasing = TRUE)]
  set.seed(seed = 1)
  enr.KEGG = GSEA(geneList = list.kegg, TERM2GENE = kegg1, eps = 0)
  
  data.plot = rbind.data.frame(slice_max(enr.KEGG@result, NES, n = 10), slice_min(enr.KEGG@result, NES, n = 10)) %>%
    mutate(Description = factor(Description, levels = Description),
           Sign = factor(sign(NES), levels = c(-1,1), labels = c("Negative", "Positive")))
  
  ggplot(data.plot, aes(x = NES, y = Description, size = -log10(`p.adjust`), col = setSize)) +
    geom_point() + scale_color_viridis() +
    facet_wrap(vars(Sign), scales = "free") +
    labs(x = "NES (Normalized enrichment score)", y = NULL,
         title = paste("Top 10 enriched KEGG terms\nfrom gene set enrichment analysis (GSEA)\non ADPDK-group",
                        cor_method, "correlations")) +
    scale_y_discrete(limits = rev) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("Output//GSEA_", cor_method, "_ADPKD-Group_Correlations.pdf"), height = 5, width = 12)
}
GSEA.correlations("Spearman")


max.data.long = max.data %>% as.data.frame %>% rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Patient", values_to = "Expression") %>%
  mutate(Group = annot.pat$Group[match(Patient, annot.pat$geo_accession)])

box.plots = map(as.list(c(data.cor$Gene[match(1:5, data.cor$Rank_Spearman)], 
                          data.cor$Gene[match(1:5, data.cor$Rank_Pearson)])), 
      function(gene){
        ggplot(filter(max.data.long, Gene == gene), aes(x = Group, y = Expression, fill = Group, group = Group)) +
          geom_boxplot() +
          labs(x = NULL, y = paste0(gene, "\nlog2 mRNA levels")) +
          theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
        })
wrap_plots(box.plots, nrow = 2) + plot_layout(guides = "collect")
ggsave("Output//Boxplots_ADPKD-Group_Correlations.pdf", height = 5, width = 15)


