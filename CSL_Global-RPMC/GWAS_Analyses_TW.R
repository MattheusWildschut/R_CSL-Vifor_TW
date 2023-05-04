setwd("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/CSL_Global-RPMC")
source("../SourceFile_TW.R")

## Downloaded from https://www.ebi.ac.uk/gwas/docs/file-downloads
data.raw = fread("Input/gwas_catalog_v1.0.2-associations_e109_r2023-03-11.tsv")
CSL = read.xlsx("Input/CSL_Global-RPMC_Overview.xlsx", sheetIndex = 1)

project = CSL$Description[2]

list.projects = map(as.list(CSL$Description), function(project){
  CSL.project = filter(CSL, Description == project)
  CSL.targets = CSL.project$Targets %>% str_split("; ") %>% map(str_remove_all, " .*") %>% unlist # %>% map(paste, collapse = " | ") %>% unlist
    # str_remove_all(" .*(?=\\;)| \\(.*\\)|\\?") %>% str_replace_all("; ", "|") %>% str_subset(".") %>% # %>% paste(collapse = "|")
    # str_split("\\|") %>% unlist
  CSL.project$Description = "GWAS associations P2RY14"
  CSL.targets = "P2RY14"
  
  if(sum(unlist(map(str_split(data.raw$MAPPED_GENE, " - |, "), function(x) any(x %in% CSL.targets)))) == 0) {
    data = matrix(nrow = 0, ncol = 0)
    } else {
      data = data.raw %>%
        # filter(str_detect(MAPPED_GENE, CSL.targets)) %>% #IFNA| ##c("P2RY14|IL3R|CSF2RB")
        filter(unlist(map(str_split(.$MAPPED_GENE, " - |, "), function(x) any(x %in% CSL.targets)))) %>%
        dplyr::mutate(MAPPED_GENE = unlist(map(str_split(MAPPED_GENE, " - |, "), function(x) paste(sort(unique(x)), collapse = " | ")))) %>%
        pivot_wider(id_cols = MAPPED_GENE, names_from = `DISEASE/TRAIT`, values_from = STUDY, values_fn = list(STUDY = length), values_fill = 0) %>%
        arrange(MAPPED_GENE) %>%
        column_to_rownames("MAPPED_GENE") %>% as.matrix
    }
  
  col_fun = colorRamp2(c(0, max(data, 1)), c("white", "red"))
  colors = structure(col_fun(0:max(data, 1)), names = as.character(0:max(data, 1)))
  draw(Heatmap(data, name = "Associations", col = colors, 
          height = unit(nrow(data)*0.5, "cm"), width = unit(ncol(data)*ifelse(ncol(data)>50, 0.3, 0.4), "cm"),
          column_title = paste0("Project: ", CSL.project$Description, "\nQueried targets: ", paste(CSL.targets, collapse = " | ")), column_title_side = "top",
          border = TRUE, rect_gp = gpar(col = "grey80"), show_row_dend = FALSE, show_column_dend = FALSE,
          column_names_gp = gpar(fontsize = 9)), padding = unit(c(6,0,0,0), "in"))
})

pdf("Output/Projects_Combined2.pdf", width = 25, height = 12)
map(list.projects, function(x){
  x
})
dev.off()

## Trial OpenTargets API --------------------------------------------------------------------------------
library(httr)
# Downloaded from https://www.genenames.org/download/custom/, selected Ensembl ID
HGNC = fread("../Other/Input/HGNC_GeneNames_04-05-23.csv")

OT_query = function(gene_id){
  ensembl_id = na.omit(HGNC$`Ensembl ID(supplied by Ensembl)`[HGNC$`Approved symbol` == gene_id])
  query_string = "
  query target($ensemblId: String!){
    target(ensemblId: $ensemblId){
      id
      approvedSymbol
      biotype
      associatedDiseases(page: { index: 0, size: 10000 }){
        count
        rows{
          disease{
            name
            therapeuticAreas{
              name
            }
            parents{
              name
              parents{
                name
                parents{
                  name
                  parents{
                    name
                    parents{
                      name
                      parents{
                        name
                        parents{
                          name
                          parents{
                            name
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          score
          datatypeScores{
            id
            score
          }
        }
      }
    }
  }
"
  
  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  variables <- list("ensemblId" = ensembl_id)
  post_body <- list(query = query_string, variables = variables)
  r <- POST(url=base_url, body=post_body, encode='json')
  
  data = httr::content(r)$data
  # x = data$target$associatedDiseases$rows[[19]]
  data3 = map_dfr(data$target$associatedDiseases$rows, function(x){
    c(setNames(c(x[[1]]$name, paste(unlist(x[[1]]$therapeuticAreas), collapse = " | "), paste(unlist(x[[1]]$parents), collapse = " | "), x[[2]]), 
               c("Disease_name", "Therapeutic_area", "Disease_parents", "Disease_score")), 
      x[[3]] %>% map_dfr(unlist) %>% pivot_wider(names_from = "id", values_from = "score"))
  }) %>% mutate(Kidney = factor(replace_na(str_extract(Disease_parents, "kidney disease|renal system measurement"), "Other"),
                                levels = c("kidney disease", "renal system measurement", "Other"), labels = c("Kidney disease", "Renal measurement", "Other"))) %>%
    group_by(Kidney) %>% slice_max(Disease_score, n = 25) %>% ungroup# %>% dplyr::select(-Kidney)
  
  data4 = data3 %>%
    dplyr::mutate(across(.cols = !Disease_name & !Therapeutic_area & !Disease_parents & !Kidney, .fns = as.numeric)) %>%
    pivot_longer(cols = -Disease_name & -Therapeutic_area & -Disease_parents & -Kidney, names_to = "Scores", values_to = "Score")
  
  areas = unique(as.vector(str_split(data4$Therapeutic_area, " \\| ", simplify = TRUE)))
  areas2 = map_dfc(as.list(areas[areas != ""]), function(x) setNames(data.frame(str_detect(data4$Therapeutic_area, x)), x)) %>%
    mutate("Disease_name" = data4$Disease_name) %>%
    pivot_longer(cols = -Disease_name) %>%
    `colnames<-`(c("Disease_name", "Area", "InArea"))
  
  diseaseTypes = c("Disease_score", "genetic_association", "somatic_mutation", "known_drug", "affected_pathway", "literature", "rna_expression", "animal_model")
  diseaseTypeLabels = c("Overall association score", "Genetic associations", "Somatic mutations", "Drugs", "Pathways & systems biology",
                        "Text mining", "RNA expression", "Animal models")
  
  data5 = data4 %>% left_join(areas2, by = "Disease_name", multiple = "all") %>%
    mutate(Disease_name = factor(Disease_name, levels = data3$Disease_name),
           # Urinary = factor(str_detect(Therapeutic_area, "urinary system disease"), levels = c(TRUE, FALSE), labels = c("Urinary", "Other")),
           Scores = factor(Scores, levels = diseaseTypes, labels = diseaseTypeLabels)) %>%
    complete(Scores, nesting(Disease_name, Therapeutic_area, Area, Kidney, InArea))
  
  plot1 = ggplot(data5, aes(x = Scores, y = Disease_name, fill = Score)) +
    geom_tile(col = "grey90") +
    scale_fill_viridis(na.value = "white") + scale_y_discrete(limits = rev) + scale_x_discrete(drop = FALSE) +
    facet_grid(rows = vars(Kidney), scales = "free", space = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = element_text(hjust = 0.5), 
          strip.background = element_blank(), strip.text.y = element_blank()) + 
    labs(x = NULL, y = NULL, title = paste("Association scores", gene_id))
  
  plot2 = ggplot(data5, aes(x = Area, y = Disease_name, fill = InArea)) +
    geom_tile(col = "grey90") +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black")) + scale_y_discrete(limits = rev) + 
    facet_grid(rows = vars(Kidney), scales = "free", space = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = element_text(hjust = 0.5),
          strip.text.y = element_text(angle = 0), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    labs(x = NULL, y = NULL, title = "Disease areas") + guides(fill = "none")
  
  plot1 + plot2 + plot_layout(guides = "collect", widths = c(1, 1.5))
}
OT_query("PADI4")
ggsave("Output/OpenTargets-Query_PADI4.pdf", width = 11, height = 10)
