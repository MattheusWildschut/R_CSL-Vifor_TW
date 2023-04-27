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
  # select(-ENTREZID) %>% 
  column_to_rownames("Gene") %>% as.matrix

library("Mfuzz")
# library("marray")

data = new("ExpressionSet", exprs = max.data)
data.s = standardise(data)
m1 = mestimate(data.s)
m1

# B) Estimate c with Dmin
# Dmin.plot = Dmin(data.s, m1, crange = seq(2,10,1), repeats = 50, visu = TRUE)
# pdf("Output//Clusters_Dmin.pdf")
# plot(Dmin.plot)
# dev.off()

# C) Estimate c with cselection
# cselection(data.s, m=m1, crange = seq(2,30,1),visu=TRUE)

#### Define Number of Clusters (c), produce and export cluster plot
c <- 6
set.seed(seed = 1)
cl <- mfuzz(data.s,c,m=m1)
mfuzz.plot2(data.s, cl=cl, mfrow = c(2,4), time.labels = colnames(raw.data),
            min.mem = 0.4, xlab = "", ylab = "",
            centre = T, centre.col = "black", centre.lwd = 3, cex.lab = 0.1, cex.axis = 1, x11 = F, srt = 90, crt = 90)

data.plot = data.frame("Gene" = rownames(data.s), "Cluster" = cl$cluster, "Membership" = rowMax(cl$membership)) %>%
  cbind.data.frame(data.s@assayData$exprs) # %>%
  # mutate(ENTREZID = str_remove(annot$ENTREZ_GENE_ID[match(Gene, annot$`Gene Symbol`)], " ///.*"))
data.plot.long = data.plot %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "Patient", values_to = "Expression") %>%
  mutate("Group" = annot.pat$Group[match(Patient, annot.pat$geo_accession)],
         "Patient2" = paste(Group, Patient))%>%
  group_by(Cluster) %>% dplyr::mutate(Cluster = paste0("Cluster ", Cluster, " (", n_distinct(Gene), " genes)"))# %>%
  # filter(Membership > 0.95)
ggplot(data.plot.long, aes(x = Patient2, y = Expression, col = Membership, group = Gene)) + #, alpha = Membership
  geom_vline(xintercept = c(3.5, 8.5, 13.5, 18.5), linetype = 2) +
  geom_line(alpha = 0.5) +
  stat_summary(fun = median, geom = "line", mapping = aes(group = -1), linewidth = 1, col = "black") + 
  scale_color_viridis(option = "B") +  #
  scale_x_discrete(limits = rev) +
  facet_wrap(. ~ Cluster, nrow = 2) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab(NULL)
# ggsave("Output\\MFuzzPlot_6cluster_Newest.pdf", width = 12, height = 8)

library(limma)
kegg2 = getKEGGPathwayNames(species="hsa")
kegg1 = getGeneKEGGLinks(species="hsa") %>%
  mutate(Gene = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, GeneID, column = "SYMBOL", keytype = "ENTREZID"),
         PathwayID = str_remove(PathwayID, "path:"),
         Description = str_remove(kegg2$Description[match(PathwayID, kegg2$PathwayID)], " - Homo sapiens \\(human\\)")) %>%
  select(Description, Gene)

kegg.data = map(as.list(1:6), function(cluster){
  list.kegg = sort(cl$membership[,cluster], decreasing = TRUE)
  set.seed(seed = 1)
  enr.KEGG = GSEA(geneList = list.kegg, scoreType = "pos", TERM2GENE = kegg1, pvalueCutoff = 1, eps = 0)
  # b = enricher(data.plot$Gene[$Cluster == 3 & data.plot$Gene %in% kegg1$Gene], TERM2GENE = kegg1, universe = data.plot$Gene)
  data.plot.cluster = slice_max(enr.KEGG@result, -log10(`p.adjust`)*sign(NES), n = 10, with_ties = FALSE) %>%
    enr.KEGG@result %>%
    mutate(Description = factor(paste(cluster, Description, sep = "_"), levels = paste(cluster, Description, sep = "_")),
           Significant = `p.adjust` < 0.05,
           Cluster = cluster)
})
data.plot = purrr::reduce(kegg.data, rbind) %>%
  mutate(Cluster = factor(Cluster, levels = 1:6, labels = paste("Cluster", 1:6)))
kegg.plots = ggplot(data.plot, aes(x = -log10(`p.adjust`), y = Description, size = NES, col = setSize, alpha = Significant)) +
  geom_point() + scale_color_viridis() +
  facet_wrap(.~Cluster, scales = "free") +
  geom_vline(xintercept = -log10(0.05), linetype = 2) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
  labs(x = "-log10(adjusted p-value)", y = NULL, #NES (Normalized enrichment score)
       title = paste0("Top 10 enriched KEGG terms\nfrom gene set enrichment analysis (GSEA)\nPer cluster (6x cluster analysis)")) +
  scale_y_discrete(limits = rev) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
ggsave("Output//Mfuzz-6Clusters_GSEA.pdf", kegg.plots, width = 17, height = 6)
