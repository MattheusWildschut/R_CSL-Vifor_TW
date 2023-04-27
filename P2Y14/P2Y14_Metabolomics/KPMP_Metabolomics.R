setwd("C:\\Users\\mwildschut\\OneDrive - Vifor Pharma AG\\Documents\\R\\Projects\\R_CSL-Vifor_TW")
source("SourceFile_TW.R")


pat.annot2 = read.csv("KPMP_Metabolomics/72ec8f42-750f-4e1b-bda9-7300d4de1365_20221208_OpenAccessClinicalData.csv")
pat.annot = read.csv("KPMP_Metabolomics/KPMP_Metabolomics_PatientAnnotations.csv", sep = ";") %>%
  `colnames<-`(str_to_title(colnames(.))) %>%
  select(-Image.type) %>%
  left_join(pat.annot2, by = join_by(Participant.id == Participant.ID)) %>%
  distinct
data = fread("KPMP_Metabolomics/metaspace_annotations_CoreMetabolome_50FDR.csv") %>%
# data = fread("KPMP_Metabolomics/metaspace_annotations_UDP_All.csv") %>%
  mutate(Patient = str_extract(datasetName, "S-[:digit:]{4}-[:digit:]{6}"),
         moleculeNames = str_replace_all(moleculeNames, ", ", "\n"),
         Database = factor(str_extract(moleculeIds, "HMDB|CHEBI|LMS"), levels = c("HMDB", "CHEBI", "LMS"),
                           labels = c("HMDB", "ChEBI", "LipidMaps"))) %>%
  mutate(moleculeNames = unlist(map(str_split(moleculeNames, "\n"), function(x) paste(sort(x), collapse = "\n")))) %>%
  mutate(moleculeNames = str_remove(moleculeNames, "\\[.*phosphonate")) %>%
  left_join(pat.annot, by = join_by(Patient == Sample.id)) %>%
  filter(!is.na(Tissue.Type))
data.plot = rbind.data.frame(data.frame("FDR" = "FDR <= 10%", filter(data, fdr<=0.1)),
                             data.frame("FDR" = "FDR <= 20%", filter(data, fdr<=0.2)),
                             data.frame("FDR" = "FDR <= 50%", filter(data, fdr<=0.5)))
b = data.plot %>% group_by(Patient, Tissue.Type) %>%
  filter(FDR == "FDR <= 10%") %>%
  dplyr::summarise(n_distinct(moleculeNames))

unique(data$moleculeNames[200])
# a = filter(data.plot, str_detect(moleculeNames, regex("uridine|UMP|UDP|UTP", ignore_case = T))) %>%
a = filter(data.plot, str_detect(moleculeNames, regex("uridine diphosphate|UDP", ignore_case = T))) %>%
  group_by(Participant.id, moleculeNames, Tissue.Type, FDR, Database) %>% dplyr::summarize(totalIntensity = mean(totalIntensity))

# a = filter(data, str_detect(moleculeNames, "glucose|Glucose|GLUCOSE"))

ggplot(a, aes(x = moleculeNames, y = totalIntensity, fill = Tissue.Type, col = Tissue.Type)) +
  geom_boxplot(position = position_dodge2(preserve = "single"), alpha = 0.25, linewidth = 0.2) +
  # guides(fill = "none") +
  facet_grid(rows = vars(FDR), scales = "free") + #, cols = vars(Database)
  # scale_y_log10() +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

## Piglet AKI model metabolomics ----------------------------------------------------------------------------------------------------
# Data downloaded from https://journals.physiology.org/doi/full/10.1152/ajprenal.00039.2022?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org

data = read.csv("KPMP_Metabolomics//Complete Metabolite Dataset_Davidson2022.csv", sep = ";") %>%
  dplyr::mutate(across(3:243, as.numeric)) %>%
  dplyr::select(1:2, contains(c("uridine", "UMP", "UDP", "UTP"))) %>%
  mutate(Group = factor(Group, labels = c("Controls", "CPB/DHCA without AKI", "CPB/DHCA with AKI")),
         Tissue = str_extract(Samples, "Serum|Kidney|Urine")) %>%
  pivot_longer(-c("Samples", "Group", "Tissue"), names_to = "Metabolite", values_to = "Intensity (Normalized)")

plots.all = ggplot(data, aes(x = Group, y = `Intensity (Normalized)`, fill = Group, group = Group)) +
  geom_boxplot(alpha = 0.5) + geom_jitter(width = 0.1, alpha = 0.5) +
  facet_grid(Tissue ~ Metabolite) +
  stat_compare_means(comparisons = list(c(1,2), c(1,3), c(2,3)), method = "t.test", method.args = list(var.equal = TRUE), label = "p.signif") +
  theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plot.udpgluc = ggplot(filter(data, Tissue == "Kidney" & Metabolite == "UDP.D.glucose"), aes(x = Group, y = `Intensity (Normalized)`, fill = Group, group = Group)) +
  geom_boxplot(alpha = 0.5) + geom_jitter(width = 0.1, alpha = 0.5) +
  facet_grid(Tissue ~ Metabolite) +
  stat_compare_means(comparisons = list(c(1,2), c(1,3), c(2,3)), method = "t.test", method.args = list(var.equal = TRUE), label = "p.signif") +
  theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plots.all + (plot_spacer() / plot.udpgluc / plot_spacer() + plot_layout(heights = c(1,1.5,1))) + plot_layout(widths = c(3,1))

             