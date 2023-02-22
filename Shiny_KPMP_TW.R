setwd("C:\\Users\\mwildschut\\OneDrive - Vifor Pharma AG\\Documents\\R\\Projects\\R_CSL-Vifor_TW")
source("SourceFile_TW.R")

# install.packages("Seurat")
# remotes::install_github("mojaveazure/seurat-disk")
# install.packages("Matrix")

library(Seurat)
library(SeuratDisk)

data = LoadH5Seurat("Input\\GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat",
                    assays = list(RNA = "data"))
data = subset(data, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
annot = fread("Input//GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_Field_Descriptions.txt.gz")
fields.pat = data.frame("Field" = c(paste0("condition", c(".long", ".l1", ".l2", ".l3")), "patient", "sex", "race")) %>%
  mutate(Description = annot$Description[match(Field, annot$Field)])
fields.cell = data.frame("Field" = c("class", paste0("subclass", c(".l1", ".l2", ".l3", ".full")), "state.l1", "state.l2")) %>%
  mutate(Description = annot$Description[match(Field, annot$Field)])

## Data from https://atlas.kpmp.org/repository/?facetTab=files&files_size=60&files_sort=%5B%7B%22field%22%3A%22
## file_name%22%2C%22order%22%3A%22asc%22%7D%5D&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in
## %22%2C%22content%22%3A%7B%22field%22%3A%22dois%22%2C%22value%22%3A%5B%2210.48698%2F92nk-e805%22%5D%7D%7D%5D%7D
# data = LoadH5Seurat("Input\\521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat")
# 
# data = subset(data, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
# fields.pat = data.frame("Field" = c("sampletype", "diseasetype", "age", "gender", "state", "tissuetype", "celltype")) #%>%
#   # mutate(Description = annot$Description[match(Field, annot$Field)])
# fields.cell = data.frame("Field" = c("ClusterClass", paste0("subclass", c(".l1", ".l2")), "celltype"))# %>%
#   # mutate(Description = annot$Description[match(Field, annot$Field)])

## UI ----------------------------------------------------------------------------------------------------------------------------------
ui = {fluidPage(
  titlePanel("CSL-Vifor_KPMP-browser_TW"),
  tabsetPanel(
    # Dotplot UI ---------
    tabPanel("Data grouping", fluid = TRUE, 
             sidebarLayout(
               sidebarPanel(
                 selectizeInput(inputId = "gene",
                                label = "Select gene of interest",
                                choices = NULL),
                 pickerInput(inputId = "pat.grouping",
                             label = "Select patient grouping", 
                             choices = fields.pat$Field,
                             choicesOpt = list(subtext = fields.pat$Description),
                             selected = "condition.l1"),
                 pickerInput(inputId = "cell.grouping",
                             label = "Select cell grouping", 
                             choices = fields.cell$Field,
                             choicesOpt = list(subtext = fields.cell$Description),
                             selected = "subclass.l2"),
                 materialSwitch(inputId = "cell.class",
                                label = "Split cell types by class",
                                value = TRUE,
                                status = "primary"),
                 uiOutput(outputId = "ui.celltype")
               ),
               mainPanel(
                 plotOutput("dotplot", height = 1000),
                 plotOutput("patient.plot")
               )
             )
    ),
  )
)}

## Server ----------------------------------------------------------------------------------------------------------------------------------
server = function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })
  updateSelectizeInput(session, "gene", choices = sort(rownames(data@assays$RNA@data)), 
                       selected = "P2RY14", options = list(maxOptions = 100000), server = TRUE)
  data.plot = reactive({
    req(input$gene, input$pat.grouping, input$cell.grouping)
    data.plot = data.frame("Expression" = data@assays$RNA@data[input$gene,], #["SLC9B2",], #
                           "Disease" = data@meta.data[,input$pat.grouping], #[,"condition.l1"], #
                           "Celltype" = data@meta.data[,input$cell.grouping], #[,"subclass.l1"], #
                           "Cellclass" = data@meta.data$class, #ClusterClass,# 
                           "Patient" = data@meta.data$patient) %>% #orig.ident) %>% # 
      group_by(Disease) %>% dplyr::mutate(Disease = paste0(Disease, " (", n_distinct(Patient), ")")) %>%
      group_by(Celltype) %>% dplyr::mutate(Celltype = paste0(Celltype, " (", length(Expression), ")"))
  })
  output$dotplot = renderPlot({
    data.plot = data.plot()
    data.sum = data.plot %>%
      group_by(Disease, Celltype, Cellclass) %>% dplyr::summarize(Percentage = sum(Expression>0)/n(), Expression = mean(Expression))

    dot.plot = ggplot(data.sum, aes(x = Disease, y = Celltype)) +
      geom_point(aes(size = Percentage, col = Expression)) +
      {if(input$cell.class) facet_grid(Cellclass ~ ., scales = "free_y", space = "free", switch = "y")} +
      # {if(input$group.pat == "Patient") facet_grid(Cellclass ~ ., scales = "free_y", space = "free", switch = "y")}
      scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(data.sum$Percentage)), guide = guide_legend(direction = "vertical")) +
      scale_color_viridis(guide = guide_colorbar(direction = "vertical")) +
      theme_bw() + theme(text = element_text(size = 18), strip.text.y.left = element_text(angle = 0), legend.position = "bottom", 
                         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    box.plot = ggplot(data.plot, aes(x = Expression, y = Celltype)) +
      # geom_violin(aes(fill = Disease)) +
      geom_boxplot(fill = "lightblue") +
      scale_x_continuous(expand = c(0,0)) + scale_fill_discrete(guide = guide_legend(direction = "vertical")) +
      {if(input$cell.class) facet_grid(Cellclass ~ ., scales = "free", space = "free", switch = "both")} +
      theme_bw() + theme(text = element_text(size = 18), legend.position = "bottom", 
                         axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), 
                         strip.background = element_blank(), strip.text.y = element_blank())
    
    dot.plot + box.plot + plot_layout(widths = c(1,3))
    # dot.plot
  })
  output$ui.celltype = renderUI({
    req(input$gene, input$pat.grouping, input$cell.grouping)
    pickerInput(inputId = "celltype",
                label = "Select single cell type of interest", 
                choices = unique(unlist(data.plot()[,"Celltype"])),
                # choicesOpt = list(subtext = fields.pat$Description),
                selected = "cDC (12)")
  })
  output$patient.plot = renderPlot({
    req(input$celltype)
    data.sum = data.plot() %>%
      group_by(Disease, Patient, Celltype) %>% dplyr::summarize(Percentage = sum(Expression>0)/n(), Expression = mean(Expression), .groups = "keep") %>%
      tibble %>% tidyr::complete(Celltype, nesting(Disease, Patient)) %>%
      filter(Celltype == input$celltype)
    ggplot(data.sum, aes(x = Patient, y = Celltype)) +
      geom_point(aes(size = Percentage, col = Expression)) +
      # {if(input$cell.class) facet_grid(Cellclass ~ ., scales = "free_y", space = "free", switch = "y")} +
      facet_grid(. ~ Disease, scales = "free", space = "free", switch = "both") +
      scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(data.sum$Percentage)), guide = guide_legend(direction = "vertical")) +
      scale_color_viridis(guide = guide_colorbar(direction = "vertical")) +
      theme_bw() + theme(text = element_text(size = 18), strip.text.y.left = element_text(angle = 0), legend.position = "bottom", 
                         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  })
  
  
  # Dotplot Server -----------
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))
