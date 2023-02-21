setwd("C:\\Users\\mwildschut\\OneDrive - Vifor Pharma AG\\Documents\\R\\Projects\\R_CSL-Vifor_TW")
source("SourceFile_TW.R")

# install.packages("Seurat")
# remotes::install_github("mojaveazure/seurat-disk")
# install.packages("Matrix")

library(Seurat)
library(SeuratDisk)

data = LoadH5Seurat("Input\\GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat",
                    assays = list(RNA = "data"))
annot = fread("Input//GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_Field_Descriptions.txt.gz")
fields.pat = data.frame("Field" = c(paste0("condition", c(".long", ".l1", ".l2", ".l3")), "patient", "sex", "race")) %>%
  mutate(Description = annot$Description[match(Field, annot$Field)])
fields.cell = data.frame("Field" = c("class", paste0("subclass", c(".l1", ".l2", ".l3", ".full")), "state.l1", "state.l2")) %>%
  mutate(Description = annot$Description[match(Field, annot$Field)])

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
                             selected = "subclass.l1"),
                 materialSwitch(inputId = "cell.class",
                                label = "Split cell types by class",
                                value = TRUE,
                                status = "primary")
               ),
               mainPanel(
                 plotOutput("dotplot", height = 800)
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
                       selected = "SLC9B2", options = list(maxOptions = 100000), server = TRUE)
  # data.plot = reactive({
  output$dotplot = renderPlot({
    req(input$gene)
    data.plot = data.frame("Expression" = data@assays$RNA@data[input$gene,], #["SLC9B2",], #
                           "Disease" = data@meta.data[,input$pat.grouping], #[,"condition.l1"], #
                           "Celltype" = data@meta.data[,input$cell.grouping], #[,"subclass.l1"], #
                           "Cellclass" = data@meta.data$class,
                           "Patient" = data@meta.data$patient) %>%
      group_by(Disease) %>% dplyr::mutate(Disease = paste0(Disease, " (", n_distinct(Patient), ")")) %>%
      group_by(Celltype) %>% dplyr::mutate(Celltype = paste0(Celltype, " (", length(Expression), ")"))
    data.sum = data.plot %>%
      group_by(Disease, Celltype, Cellclass) %>% dplyr::summarize(Expression = mean(Expression), Percentage = sum(Expression>0)/n())

    dot.plot = ggplot(data.sum, aes(x = Disease, y = Celltype)) +
      geom_point(aes(size = Percentage, col = Expression)) +
      {if(input$cell.class) facet_grid(Cellclass ~ ., scales = "free_y", space = "free", switch = "y")} +
      scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(data.sum$Percentage))) +
      scale_color_viridis() +
      theme_bw() + theme(text = element_text(size = 18), strip.text.y.left = element_text(angle = 0), legend.position = "bottom")
    
    violin.plot = ggplot(data.plot, aes(x = Expression, y = Celltype)) +
      geom_violin(aes(fill = Disease)) +
      scale_x_continuous(expand = c(0,0)) +
      {if(input$cell.class) facet_grid(Cellclass ~ ., scales = "free_y", space = "free", switch = "y")} +
      theme_bw() + theme(text = element_text(size = 18), legend.position = "bottom", axis.text.y = element_blank(), axis.title.y = element_blank(),
                         strip.background = element_blank(), strip.text.y = element_blank())
    
    dot.plot + violin.plot
  })
  # output$dotplot = renderPlot({
    
  # })
  
  
  # Dotplot Server -----------
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))
