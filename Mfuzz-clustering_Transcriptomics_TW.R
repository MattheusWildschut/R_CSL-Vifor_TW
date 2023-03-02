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
Dmin.plot = Dmin(data.s, m1, crange = seq(2,10,1), repeats = 50, visu = TRUE)
plot(Dmin.plot)

# C) Estimate c with cselection
cselection(data.s, m=m1, crange = seq(2,30,1),visu=TRUE)

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
ggsave("Output\\MFuzzPlot_6cluster_Newest.pdf", width = 12, height = 8)

library(limma)
kegg2 = getKEGGPathwayNames(species="hsa")
kegg1 = getGeneKEGGLinks(species="hsa") %>%
  mutate(Gene = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, GeneID, column = "SYMBOL", keytype = "ENTREZID"),
         PathwayID = str_remove(PathwayID, "path:"),
         Description = str_remove(kegg2$Description[match(PathwayID, kegg2$PathwayID)], " - Homo sapiens \\(human\\)")) %>%
  select(Description, Gene)

kegg.plots = map(as.list(1:6), function(cluster){
  set.seed(seed = 1)
  c = data.plot %>%
    filter(Gene %in% kegg1$Gene & Cluster == cluster) %>%
    arrange(desc(Membership))
  d = setNames(c$Membership, c$Gene)
  e = GSEA(geneList = d, scoreType = "pos", TERM2GENE = kegg1, pvalueCutoff = 0.1)
  # b = enricher(data.plot$Gene[data.plot$Cluster == 3 & data.plot$Gene %in% kegg1$Gene], TERM2GENE = kegg1, universe = data.plot$Gene)
  ggplot(e@result[1:10,], aes(x = -log10(`p.adjust`), y = reorder(Description, `p.adjust`))) +#, size = Count, col = GeneRatio)) +
    geom_point() + #scale_color_viridis() +
    # facet_wrap(vars(Sign), scales = "free") +
    labs(y = NULL, #x = "NES (Normalized enrichment score)", 
         title = paste("Top 10 enriched KEGG terms\nCluster", cluster)) + #\non ", 
    # input$gene, " (", ifelse(input$probe.sum, "Max probe", select_probe()), ") correlations")) +
    scale_y_discrete(limits = rev) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
})
wrap_plots(kegg.plots, nrow = 2) + plot_layout(guides = "collect")

gene.ID = uniprot$GeneID[match(data.plot$Gene, uniprot$`Gene Names (primary)`)]
gene.ID = str_remove(gene.ID, ";.*")
gene.cluster = gene.ID[data.plot$Cluster == 3 & data.plot$Membership > 0.95 & !is.na(gene.ID) & gene.ID != ""]
set.seed(1)
k <- kegga(gene.cluster, universe = gene.ID[!is.na(gene.ID) & gene.ID != ""],
           species="Hs", plot = TRUE, pathway.names = getKEGGPathwayNames("hsa", remove=TRUE))
top.k = topKEGG(k)
top.k$Pathway = str_remove(kegg2$Description[match(str_remove(rownames(top.k), "path:"), kegg2$PathwayID)], " - Homo sapiens \\(human\\)")

ggplot(top.k[1:10,], aes(x = -log10(P.DE), y = reorder(Pathway, P.DE), size = N, col = P.DE)) +
  geom_point() + #scale_color_viridis() +
  # facet_wrap(vars(Sign), scales = "free") +
  labs(y = NULL, #x = "NES (Normalized enrichment score)", 
       title = paste0("Top 10 enriched KEGG terms\nfrom gene set enrichment analysis (GSEA)")) + #\non ", 
  # input$gene, " (", ifelse(input$probe.sum, "Max probe", select_probe()), ") correlations")) +
  scale_y_discrete(limits = rev) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))



length(unique(kegg1$GeneID))
a = data.plot %>% group_by(Cluster) %>% slice_max(Membership, n = 1000)

z = KEGGREST::keggConv("ncbi-geneid", "hsa")
y = names(z)[match(data.plot$Gene, str_remove(z, "ncbi-geneid:"))]
y[1:10]
data.plot$Gene[1:10]

kegg2 = getKEGGPathwayNames(species="hsa")
kegg1 = getGeneKEGGLinks(species="hsa") %>%
  mutate(Gene = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, GeneID, column = "SYMBOL", keytype = "ENTREZID"),
         PathwayID = str_remove(PathwayID, "path:"),
         Description = str_remove(kegg2$Description[match(PathwayID, kegg2$PathwayID)], " - Homo sapiens \\(human\\)")) %>%
  select(Description, Gene)

data.plot$Gene[(!data.plot$Gene %in% kegg1$Gene)]
unique(kegg1$Gene)

kegg1$Description = kegg2$Description[match(kegg1$PathwayID, kegg2$PathwayID)]
kegg1 = cbind.data.frame(kegg1$Description, kegg1$GeneID)
# tab = cbind.data.frame(tab$)
# tab$Symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, tab$GeneID,
#                        column="SYMBOL", keytype="ENTREZID")
b = enricher(data.plot$Gene[data.plot$Cluster == 2], TERM2GENE = kegg1) #universe = data.plot$Gene)
dotplot(b)

ggplot(b@result[1:10,], aes(x = -log10(`p.adjust`), y = reorder(Description, `p.adjust`), size = Count, col = GeneRatio)) +
  geom_point() + #scale_color_viridis() +
  # facet_wrap(vars(Sign), scales = "free") +
  labs(y = NULL, #x = "NES (Normalized enrichment score)", 
       title = paste0("Top 10 enriched KEGG terms\nfrom gene set enrichment analysis (GSEA)")) + #\non ", 
                      # input$gene, " (", ifelse(input$probe.sum, "Max probe", select_probe()), ") correlations")) +
  scale_y_discrete(limits = rev) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))



d = sort(cl$membership[,4], decreasing = TRUE)
e = scale(d)[,1]
names(e) = str_remove(annot$ENTREZ_GENE_ID[match(names(e), annot$`Gene Symbol`)], " ///.*")
e = e[!duplicated(names(e))]

f = gseKEGG(e, keyType = "ncbi-geneid", organism = "hsa", eps = 0, pvalueCutoff = 1)
g = f@result
colnames(b@result)
