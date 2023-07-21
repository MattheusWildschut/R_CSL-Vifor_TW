setwd("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/Other")
source("../SourceFile_TW.R")

data = readxl::read_xlsx("Input/Hongbo Liu_2023_S24.xlsx", sheet = "Table S24", skip = 2) %>%
  pivot_longer(cols = contains("Specific Peak"), names_to = "CellType", values_to = "Specific") %>%
  filter(Specific == "Yes") %>%
  mutate(CellType = str_remove(CellType, "-Specific Peak"))# %>%
  # filter(`Priority Score` == 8)

write_tsv(data, "Output/Susztak_KidneyGenes_S24_All.tsv")

node_att = data.frame("node_id" = c(unique(data$`Top Targe Gene of Locus`), unique(data$CellType)),
                      "node_type" = rep(c("Gene", "CellType"), c(length(unique(data$`Top Targe Gene of Locus`)), length(unique(data$CellType)))))
write_tsv(node_att, "Output/Susztak_KidneyGenes_S24_NodeAttributes.tsv")
