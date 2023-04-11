Sys.setenv(LANG = "EN")
library(readxl)
library(xlsx)
library(stringr)
library(dplyr)
library(purrr)

## Merge excel files for different folders ------------------------------------------------------------------------
root.folder = "C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/Projects/RESOLUTE/"
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
root.folder = "C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/Projects/RESOLUTE/"
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
