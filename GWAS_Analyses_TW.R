setwd("C:\\Users\\mwildschut\\OneDrive - Vifor Pharma AG\\Documents\\R\\Projects\\R_CSL-Vifor_TW")
source("SourceFile_TW.R")

## Downloaded from https://www.ebi.ac.uk/gwas/docs/file-downloads
data.raw = fread("Input/gwas_catalog_v1.0.2-associations_e109_r2023-03-11.tsv")
CSL = read.xlsx("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/Projects/CSL_Global-RPMC/CSL_Global-RPMC_Overview.xlsx",
                sheetIndex = 1)

project = CSL$Description[3]

list.projects = map(as.list(CSL$Description), function(project){
  CSL.project = filter(CSL, Description == project)
  CSL.targets = CSL.project$Targets %>% str_replace_all("; ", "|") %>% str_remove_all(" .*|\\?") %>% str_subset(".") %>% # %>% paste(collapse = "|")
    str_split("\\|") %>% unlist
  
  if(sum(unlist(map(str_split(data.raw$MAPPED_GENE, " - |, "), function(x) any(x %in% CSL.targets)))) == 0) {
    data = matrix(nrow = 0, ncol = 0)
    } else {
      data = data.raw %>%
        # filter(str_detect(MAPPED_GENE, CSL.targets)) %>% #IFNA| ##c("P2RY14|IL3R|CSF2RB")
        filter(unlist(map(str_split(.$MAPPED_GENE, " - |, "), function(x) any(x %in% CSL.targets)))) %>%
        dplyr::mutate(MAPPED_GENE = unlist(map(str_split(MAPPED_GENE, " - |, "), function(x) paste(sort(x), collapse = " | ")))) %>%
        pivot_wider(id_cols = MAPPED_GENE, names_from = `DISEASE/TRAIT`, values_from = STUDY, values_fn = list(STUDY = length), values_fill = 0) %>%
        arrange(MAPPED_GENE) %>%
        column_to_rownames("MAPPED_GENE") %>% as.matrix
    }
  
  col_fun = colorRamp2(c(0, max(data, 1)), c("white", "red"))
  colors = structure(col_fun(0:max(data, 1)), names = as.character(0:max(data, 1)))
  draw(Heatmap(data, name = "Associations", col = colors, 
          height = unit(nrow(data)*0.5, "cm"), width = unit(ncol(data)*ifelse(ncol(data)>50, 0.3, 0.4), "cm"),
          column_title = paste0("Project: ", CSL.project$Description, "\nQueried targets: ", paste(CSL.targets, collapse = " | ")), column_title_side = "top",
          border = TRUE, rect_gp = gpar(col = "grey80"), 
          column_names_gp = gpar(fontsize = 9)), padding = unit(c(6,0,0,0), "in"))
})

pdf("CSL_Global-RPMC/Output/Projects_Combined3.pdf", width = 25, height = 12)
map(list.projects, function(x){
  x
})
dev.off()
