setwd("C:\\Users\\mwildschut\\OneDrive - Vifor Pharma AG\\Documents\\R\\Projects\\R_CSL-Vifor_TW")
source("SourceFile_TW.R")

library(GEOquery)
dataset = "GSE7869" #"GSE35831" #
gse = getGEO(dataset, GSEMatrix = TRUE)[[1]]

data = gse@assayData$exprs
annot = gse@featureData@data
annot.pat = gse@phenoData@data %>%
  mutate(Group = factor(case_when(str_detect(title, "normal|Normal") ~ "Normal",
                                  str_detect(title, "minimally") ~ "MCD",
                                  str_detect(title, "small") ~ "Cysts-Small",
                                  str_detect(title, "medium") ~ "Cysts-Medium",
                                  str_detect(title, "large") ~ "Cysts-Large",
                                  str_detect(title, "ADPKD") ~ "ADPKD"),
                        levels = c("Normal", "MCD", "Cysts-Small", "Cysts-Medium", "Cysts-Large", "ADPKD")))
if(dataset == "GSE35831"){
  annot$`Gene Symbol` = str_split_fixed(annot$gene_assignment, " // ", 4)[,2]
  annot$`Gene Title` = str_split_fixed(annot$gene_assignment, " // ", 4)[,3]
}
IDs = as.character(annot$ID[annot$`Gene Symbol` == "SLC9B2"])
probe = 1
cors = as.data.frame(cor(t(data), data[IDs[probe],])) %>%
  `colnames<-`("Correlation") %>%
  mutate(ID = rownames(.),
         Gene = annot$`Gene Symbol`[match(annot$ID, ID)],
         Description = annot$`Gene Title`[match(annot$ID, ID)]) %>%
  filter(!rownames(.) == IDs[probe]) %>%
  arrange(-Correlation)

cor_function = function(direction, i){
  if(direction == "neg"){
    cors.sel = cors %>%
      arrange(Correlation)
  } else {cors.sel = cors}
  idx = cors.sel[i,]
  data.cor = cbind.data.frame(data[IDs[probe],], data[idx$ID,]) %>%
    `colnames<-`(c("SLC9B2", "Gene")) %>%
    mutate(Group = annot.pat$Group[match(annot.pat$geo_accession, rownames(.))])
  scatter.plot = ggplot(data.cor, aes(x = SLC9B2, y = Gene)) +
    geom_point(aes(col = Group)) +
    geom_smooth(method = "lm", formula = y ~ x) + stat_cor() +
    labs(x = paste0("SLC9B2 (", IDs[probe], ")"), y = paste0(idx$Gene, " (", idx$ID, ")")) +
    theme_bw()
  box.plot = ggplot(data.cor, aes(x = Group, y = Gene)) +
    geom_boxplot(aes(fill = Group)) +
    labs(x = NULL, y = paste0(idx$Gene, " (", idx$ID, ")")) +
    theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  return(list("scatter.plot" = scatter.plot, "box.plot" = box.plot))
}
scatter.list = map(1:5, function(x) cor_function("pos", x)$scatter.plot)
boxplot.list = map(1:5, function(x) cor_function("pos", x)$box.plot)
wrap_plots(scatter.list, nrow = 1) / wrap_plots(boxplot.list, nrow = 1) + plot_layout(guides = "collect")

cors.sum = cors %>%
  mutate(ENTREZID = str_remove(annot$ENTREZ_GENE_ID[match(ID, annot$ID)], " ///.*")) %>%
  filter(ENTREZID != "") %>%
  group_by(ENTREZID) %>%
  dplyr::summarise(Correlation = mean(Correlation)) %>%
  arrange(-Correlation)

# enr.pos = enrichKEGG(gene = slice_max(cors.sum, Correlation, n = 250)$ENTREZID, universe = cors.sum$ENTREZID)
# enr.neg = enrichKEGG(gene = slice_min(cors.sum, Correlation, n = 250)$ENTREZID, universe = cors.sum$ENTREZID)
# dotplot(enr.pos) + dotplot(enr.neg)

cors.list = setNames(cors.sum$Correlation, cors.sum$ENTREZID)
set.seed(0)
enr.KEGG = gseKEGG(cors.list, keyType = "ncbi-geneid", organism = "hsa", eps = 0)
# dotplot(enr.KEGG, x = "p.adjust", showCategory=10, split=".sign") + facet_wrap(.~.sign, scales = "free_y")
direction = "neg"
data.plot = rbind.data.frame(slice_max(enr.KEGG@result, NES, n = 10), slice_min(enr.KEGG@result, NES, n = 10)) %>%
  mutate(Description = factor(Description, levels = Description),
         Sign = factor(sign(NES), levels = c(-1,1), labels = c("Negative", "Positive")))

ggplot(data.plot, aes(x = NES, y = Description, size = -log10(`p.adjust`), col = setSize)) +
  geom_point() + scale_color_viridis() +
  facet_wrap(vars(Sign), scales = "free") +
  labs(x = "NES (Normalized enrichment score)", y = NULL, 
       title = paste("Top 10 enriched KEGG terms\nfrom gene set enrichment analysis (GSEA)\non NHA2/SLC9B2 probe", IDs[probe], "correlations")) +
  scale_y_discrete(limits = rev) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("Output\\DotPlot_GSEA_NHA2-Correlations_Song2009_Probe", probe, ".pdf"), width = 12, height = 6, units = "in")

gseaplot(enr.KEGG, geneSetID = 3, by = "runningScore", title = edox$Description[3])
a = enr.KEGG@result
