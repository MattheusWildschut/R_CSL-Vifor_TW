## Merge excel files for different folders ------------------------------------------------------------------------
Sys.setenv(LANG = "EN")
library(readxl)
library(xlsx)
library(stringr)
library(purrr)
library(dplyr)

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
