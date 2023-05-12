data = data.snRNA
data@meta.data = data@meta.data %>%
  mutate(Celltype = subclass.l1,
         Disease = diseasetype)
celltypes = unique(data@meta.data$Celltype)

data.pct.pat = Percent_Expressing(data, features = "AHSG",
                                  group_by = "Patient", 
                                  split_by = "Celltype") %>%
  `colnames<-`(str_replace(str_remove(colnames(.), "^X"), "\\.", "/")) %>%
  t %>% as.data.frame %>%
  mutate(Patient = str_split_fixed(rownames(.), "_", 2)[,1],
         Celltype = str_split_fixed(rownames(.), "_", 2)[,2]) %>%
  mutate(Percentage = AHSG/100)

data.counts.pat = data@meta.data %>%
  dplyr::count(Patient, Celltype, .drop = FALSE, name = "CellCounts") %>%
  left_join(dplyr::count(data@meta.data, Patient, name = "AllCells"), by = "Patient") %>%
  left_join(data@meta.data %>% group_by(Disease) %>% dplyr::mutate("Patients" = n_distinct(Patient)) %>%
              dplyr::select(Patient, Patients) %>% distinct, by = "Patient") %>%
  full_join(data.frame("Celltype" = rep(celltypes, each = length(unique(data@meta.data$Patient))),
                       "Patient" = unique(data@meta.data$Patient)), by = c("Patient", "Celltype"), multiple = "all")

data.avg.pat = AverageExpression(subset(data, features = "AHSG"),
                                 group.by = c("Patient", "Celltype"))[[1]]

data.patient = data.avg.pat %>%
  t %>% as.data.frame %>% `colnames<-`("Expression") %>%
  mutate(Patient = str_split_fixed(rownames(.), "_", 2)[,1],
         Celltype = str_split_fixed(rownames(.), "_", 2)[,2]) %>%
  full_join(data.pct.pat, by = c("Celltype", "Patient")) %>%
  full_join(data.counts.pat, by = c("Celltype", "Patient")) %>%
  mutate(Patient = paste0(Patient, " (", AllCells, ")"),
         Disease = paste0(Disease, " (", Patients, ")"))

ggplot(data.patient, aes(x = reorder(Patient, -CellCounts), y = Celltype)) +
  geom_point(aes(size = Percentage, col = Expression)) +
  facet_grid(. ~Disease, scales = "free", space = "free", switch = "x") +
  scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(c(0, data.patient$Percentage), na.rm = TRUE)), guide = guide_legend(direction = "vertical")) +
  scale_color_viridis(guide = guide_colorbar(direction = "vertical")) +
  theme_bw() + theme(text = element_text(size = 18), strip.text.y.left = element_text(angle = 0), legend.position = "bottom",
                     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("Patient")
ggsave("Output/KPMP-AHSG_snRNAseq_Celltypes.pdf", width = 9, height = 7)

pat.annot3 = pat.annot2 %>% select(Participant.ID, Baseline.eGFR..ml.min.1.73m2...Binned.) %>%
  `colnames<-`(c("Patient", "eGFR")) %>%
  mutate(Patient = str_remove(Patient, "-"),
         eGFR = as.numeric(str_remove(eGFR, "-.*")))
data.patient2 = data.patient %>%
  mutate(Patient = str_extract(str_remove(Patient, " \\(.*"), "[:digit:]+")) %>%
  left_join(pat.annot3, by = "Patient") %>%
  filter(!is.na(eGFR))

plot.pct = ggplot(data.patient2, aes(x = Percentage, y = eGFR, col = Disease)) +
  geom_point() +
  facet_wrap(.~Celltype, scales = "free") + theme_bw()
plot.expr = ggplot(data.patient2, aes(x = Expression, y = eGFR, col = Disease)) +
  geom_point() +
  facet_wrap(.~Celltype, scales = "free") + theme_bw()
plot.pct + plot.expr + plot_layout(guides = "collect")
ggsave("Output/KPMP-AHSG_snRNAseq_eGFR.pdf", width = 15, height = 7)
