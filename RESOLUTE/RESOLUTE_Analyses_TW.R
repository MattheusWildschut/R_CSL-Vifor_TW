Sys.setenv(LANG = "EN")
library(readxl)
library(xlsx)
library(stringr)
library(dplyr)
library(purrr)

## Merge excel files for different folders ------------------------------------------------------------------------
root.folder = "C:/Users/mwildschut/OneDrive - CSL Behring/Documents/Projects/RESOLUTE/"
folders = list("extrapolated_concentrations_results", "intensity_results", "concentration_results", "experiments_with_only_conc_analyzed")

walk(folders, function(folder){
  print(paste("Summarizing folder", folder))
  setwd(paste0(root.folder, "Ionomics/DATA/summary/"))
  setwd(folder)
  files = list.files(pattern = ".xlsx")
  file.names = str_extract(files, "IONEX[:digit:]*")
  
  tab.names = excel_sheets(files[1])
  wb = createWorkbook()
  map(as.list(tab.names), function(tab){
    print(paste("Adding tab", tab))
    data.list = map2(as.list(files), as.list(file.names), function(x,y){
      if(tab %in% excel_sheets(x)){
        data.frame("resoluteID" = y,
                   read.xlsx(x, sheetName = tab) %>% filter(rowSums(is.na(.)) != ncol(.))) %>%
          `colnames<-`(str_replace(colnames(.), "Tretment|Treatement", "Treatment")) %>%
          `colnames<-`(ifelse(str_detect(colnames(.), tab), colnames(.), paste0(tab, ".", colnames(.))))
      }
    })
    data.merged = reduce(data.list, bind_rows)
    addDataFrame(data.merged, sheet = createSheet(wb, tab), row.names = FALSE)
  })
  setwd(paste0(root.folder, "Ionomics/DATA/summary/Merged_TW"))
  saveWorkbook(wb, paste0("Merged_", folder, "_TW.xlsx"))
})

## Differential results merging ------------------------------------------------------------------------------------------
root.folder = "C:/Users/mwildschut/OneDrive - CSL Behring/Documents/Projects/RESOLUTE/"
folders = list("extrapolated_concentrations_results", "experiments_with_only_conc_analyzed", "concentration_results",  "intensity_results")

setwd(paste0(root.folder, "Ionomics/DATA/summary/Merged_TW/"))
tests.list = map(folders, function(folder){
  print(paste("Summarizing", folder))
  Metadata = read.xlsx(paste0("Merged_", folder, "_TW.xlsx"), sheetName = "Metadata")
  tests = read.xlsx(paste0("Merged_", folder, "_TW.xlsx"), sheetName = "tests") %>%
    `colnames<-`(str_remove(colnames(.), "tests.")) %>%
    mutate(dataset_file_name = paste(str_extract(experiment, "^[:digit:]*"), "Resolute", resoluteID, sep = "_"),
           id = paste(ion, testName, dataset_file_name, sep = "_"),
           dashboard_value = str_replace(folder, "_results", "_value"),
           ExpName = Metadata$Metadata.ExpName[match(resoluteID, Metadata$Metadata.resoluteID)]) %>%
    select(id, dashboard_value, dataset_file_name, ExpName, testName, ion, lfc, pval, padj, padj.signif) %>%
    `colnames<-`(c("id", "dashboard_value", "dataset_file_name", "ExpName", "SLCShort", "ion", "log2FC", "pvalue", "padj", "padj.signif"))
})
tests.merged = purrr::reduce(tests.list, bind_rows) %>%
  arrange(id)
tests.filtered = tests.merged %>% 
  mutate(dashboard_value = factor(dashboard_value, levels = str_replace(unlist(folders), "_results", "_value"))) %>%
  group_by(id) %>% slice_min(dashboard_value, n = 1) %>% as.data.frame
table(tests.filtered$dashboard_value)

wb = createWorkbook()
addDataFrame(tests.filtered, sheet = createSheet(wb, "Filtered"), row.names = FALSE)
addDataFrame(tests.merged, sheet = createSheet(wb, "All"), row.names = FALSE)
saveWorkbook(wb, "Data_to_visualize_on_dashboard_TW.xlsx")

## Make Excel of selected SLCs for Alvaro ------------------------------------------------------------------------------
setwd("C:/Users/mwildschut/OneDrive - CSL Behring/Documents/Projects/RESOLUTE/Ionomics/DATA/summary/Merged_TW/")
SLCs.sel = c("ANKH", "NIPAL2", "SLC39A11", "MFSD3", "SLC16A6", "SLC16A4", "SLC35F1", "SLC16A13", "SLC45A3", "SLCO5A1", "MFSD5", "SLC16A5", "SLC45A4")
data.raw = read.xlsx("Data_to_visualize_on_dashboard_TW.xlsx", sheetName = "Filtered") 
data.sel = data.raw %>%
  filter(SLCShort %in% SLCs.sel) %>%
  filter(padj < 0.05) %>%
  mutate(SLCShort = factor(SLCShort, levels = SLCs.sel)) %>%
  arrange(SLCShort) %>%
  merge(data.frame("SLCShort" = SLCs.sel), by = "SLCShort", all = TRUE)
wb = createWorkbook()
addDataFrame(data.sel, sheet = createSheet(wb, "Selected_SLCs"), row.names = FALSE)
saveWorkbook(wb, "Selected_SLCs_ForAlvaro_TW.xlsx")

## Make plot for presentation
setwd("C:/Users/mwildschut/OneDrive - CSL Behring/Documents/Projects/RESOLUTE/Ionomics/DATA/summary/Merged_TW/")

data1 = read.xlsx("Merged_extrapolated_concentrations_results_TW.xlsx", sheetName = "Data")
data2 = read.xlsx("Merged_extrapolated_concentrations_results_TW.xlsx", sheetName = "tests")

data.plot2 = data2 %>% 
  filter(str_detect(tests.testName, "SLC39") & !str_detect(tests.testName, "SLC39A6") &
           !(tests.testName == "SLC39A14" & tests.experiment %in% c("202102_Resolute_IONEX0015", "211027_Resolute_IONEX0029"))) %>%
  mutate(tests.ion = str_remove(tests.ion, "\\(.*"), SLC = as.numeric(str_extract(tests.testName, "[:digit:]+$"))) %>%
  arrange(SLC)

library(ggplot2)
library(viridis)
library(scales)
library(RColorBrewer)
ggplot(data.plot2, aes(x = tests.ion, y = reorder(tests.testName, -SLC), size = -log10(tests.padj), fill = tests.lfc)) +
  geom_point(shape = 21, col = ifelse(data.plot2$tests.padj < 0.05, "black", "grey")) +
  colorspace::scale_fill_continuous_divergingx(palette = "RdBu", mid = 0) +
  scale_size_continuous(breaks = -log10(c(0.99, 0.5, 0.05, 0.0001, 0.000001)), labels = c(1, 0.5, 0.05, 0.0001, 0.000001)) +
  labs(x = "Measured ion", y = "SLC39 family", fill = "Log2(fold change)", size = "Adjusted p-value") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                     panel.grid = element_blank())

data.plot1 = data1 %>%
  filter(Data.SLCShort == "SLC39A14" & Data.ion == "Fe56(MR)" & Data.ExpName == "IA_36") %>%
  mutate(Condition = factor(Data.Treatment, levels = c("Vehicle", "+dox (1 µg/ml)"), labels = c("Control (DOX-)", "Induced (DOX+)")))

ggplot(data.plot1, aes(x = Condition, y = Data.intensity, col = Condition, fill = Condition)) +
  geom_boxplot(alpha = 0.25) + geom_jitter(width = 0.2, shape = 21, col = "white") +
  labs(x = NULL, y = "Fe56 ion concentration (µg/L)\nnormalized to protein (mg/mL)") +
  theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     panel.grid = element_blank())

data3 = read.xlsx("Data_to_visualize_on_dashboard_TW.xlsx", sheetName = "All", header = TRUE)
unique(str_remove(data3$SLCShort, "mito|_.*"))

a = c(read.xlsx("Merged_concentration_results_TW.xlsx", sheetName = "Metadata") %>% dplyr::select(Metadata.SLCShort),
      read.xlsx("Merged_extrapolated_concentrations_results_TW.xlsx", sheetName = "Metadata") %>% dplyr::select(Metadata.SLCShort),
      read.xlsx("Merged_experiments_with_only_conc_analyzed_TW.xlsx", sheetName = "Metadata") %>% dplyr::select(Metadata.SLCShort),
      read.xlsx("Merged_intensity_results_TW.xlsx", sheetName = "Metadata") %>% dplyr::select(Metadata.SLCShort))
sort(unique(str_remove(unlist(a), "mito|_.*")))

data.annot = read.xlsx("C:/Users/mwildschut/OneDrive - CSL Behring/Documents/R/Projects/R_CSL-Vifor_TW/RESOLUTE/Input/Copy of ionomics vifor pri.xlsx",
                       sheetIndex = 1, startRow = 8)

data.merged = merge(data.frame(Measured = "Measured", SLC.Name = unique(str_remove(data3$SLCShort, "mito|_.*"))),
                    data.annot, by = "SLC.Name", all = TRUE)
write.csv(data.merged, "C:/Users/mwildschut/OneDrive - CSL Behring/Documents/R/Projects/R_CSL-Vifor_TW/RESOLUTE/Input/Annotated_SLCs.csv")
