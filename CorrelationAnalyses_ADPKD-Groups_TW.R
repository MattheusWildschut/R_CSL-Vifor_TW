## Data loading
setwd("C:\\Users\\mwildschut\\OneDrive - Vifor Pharma AG\\Documents\\R\\Projects\\R_CSL-Vifor_TW")
source("SourceFile_TW.R")

dataset = "GSE7869"
gse = getGEO(dataset, GSEMatrix = TRUE)[[1]]

uniprot = fread("Input\\Uniprot_HumanProteome_28-02-23.tsv")
annot = gse@featureData@data
raw.data = gse@assayData$exprs[annot$`Gene Symbol` != "" & !str_detect(annot$`Gene Symbol`, "-| ///|\\."),]
# & annot$`Gene Symbol` %in% uniprot$`Gene Names (primary)`,]
annot = annot[annot$`Gene Symbol` != "" & !str_detect(annot$`Gene Symbol`, "-| ///|\\."),]
# & annot$`Gene Symbol` %in% uniprot$`Gene Names (primary)`,]
annot.pat = gse@phenoData@data %>%
  mutate(Group = factor(case_when(str_detect(title, "normal|Normal") ~ "Normal",
                                  str_detect(title, "minimally") ~ "MCD",
                                  str_detect(title, "small") ~ "Cysts-Small",
                                  str_detect(title, "medium") ~ "Cysts-Medium",
                                  str_detect(title, "large") ~ "Cysts-Large"),
                        levels = c("Normal", "MCD", "Cysts-Small", "Cysts-Medium", "Cysts-Large")))

annot.genes = data.frame("Gene" = sort(unique(annot$`Gene Symbol`)), 
                         "Description" = annot$`Gene Title`[match(sort(unique(annot$`Gene Symbol`)), annot$`Gene Symbol`)])

max.data = as.data.frame(raw.data) %>%
  mutate(Gene = annot$`Gene Symbol`) %>%
  # ENTREZID = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, Gene, keytype="SYMBOL", column = "ENTREZID")) %>%
  # filter(!is.na(ENTREZID)) %>%
  group_by(Gene) %>%
  dplyr::summarise(across(everything(), max)) %>%
  column_to_rownames("Gene") %>% as.matrix
  
max.data.long = max.data %>%
  pivot_longer(-Gene, names_to = "Patient", values_to = "Expression") %>%
  mutate(Group = as.numeric(annot.pat$Group[match(Patient, annot.pat$geo_accession)]))
  # select(-ENTREZID) %>% 
  # column_to_rownames("Gene") %>% as.matrix

a = foreach(i = unique(max.data.long$Gene), .combine = rbind) %dopar% {
  res.lm = lm(data = filter(max.data.long, Gene == i), formula = Group ~ Expression)
  sum.lm = summary(res.lm)
  res.aov = Anova(res.lm)
  c(i, -log10(res.aov$`Pr(>F)`[-nrow(res.aov)])*ifelse(sum.lm[[4]][,"Estimate"][-1] > 0, 1, -1))
}

b = cor(t(max.data), as.numeric(annot.pat$Group), method = "pearson") %>% 
  as.data.frame %>% dplyr::rename(Correlation = V1) %>% rownames_to_column("Gene") %>% arrange(desc(Correlation))
c = setNames(b$Correlation, b$Gene)

library(limma)
kegg2 = getKEGGPathwayNames(species="hsa")
kegg1 = getGeneKEGGLinks(species="hsa") %>%
  mutate(Gene = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, GeneID, column = "SYMBOL", keytype = "ENTREZID"),
         PathwayID = str_remove(PathwayID, "path:"),
         Description = str_remove(kegg2$Description[match(PathwayID, kegg2$PathwayID)], " - Homo sapiens \\(human\\)")) %>%
  select(Description, Gene)

set.seed(seed = 1)
enr.KEGG = GSEA(geneList = c, TERM2GENE = kegg1, eps = 0)

data.plot = rbind.data.frame(slice_max(enr.KEGG@result, NES, n = 10), slice_min(enr.KEGG@result, NES, n = 10)) %>%
  mutate(Description = factor(Description, levels = Description),
         Sign = factor(sign(NES), levels = c(-1,1), labels = c("Negative", "Positive")))

ggplot(data.plot, aes(x = NES, y = Description, size = -log10(`p.adjust`), col = setSize)) +
  geom_point() + scale_color_viridis() +
  facet_wrap(vars(Sign), scales = "free") +
  labs(x = "NES (Normalized enrichment score)", y = NULL,
       title = paste0("Top 10 enriched KEGG terms\nfrom gene set enrichment analysis (GSEA)\non ADPDK-group correlations")) +#, 
                      # input$gene, " (", ifelse(input$probe.sum, "Max probe", select_probe()), ") correlations")) +
  scale_y_discrete(limits = rev) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
ggsave("Output//GSEA_Pearson_ADPKD-Group_Correlations.pdf", height = 5, width = 12)
