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
                 uiOutput(outputId = "ui.celltype"),
                 downloadBttn('download.xlsx', list(icon("file-excel"), "Download data as Excel", icon("file-excel")),
                              color = "success", size = "sm")
               ),
               mainPanel(
                 plotOutput("plot.celltype", height = 1000),
                 uiOutput("ui.plot.patient")
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
  output$ui.celltype = renderUI({
    celltypes = unique(sort(as.character(unlist(data@meta.data[,input$cell.grouping]))))
    pickerInput(inputId = "celltype",
                label = "Select cell types of interest for single patient assesment",
                choices = celltypes,
                multiple = TRUE,
                selected = NULL,
                options = list(`actions-box` = TRUE))
  })
  data.celltype = reactive({
    req(input$gene, input$pat.grouping, input$cell.grouping)
    data.plot = data.frame("Expression" = data@assays$RNA@data[input$gene,], #["SLC9B2",], #
                           "Disease" = data@meta.data[,input$pat.grouping], #[,"condition.l1"], #
                           "Celltype" = data@meta.data[,input$cell.grouping], #[,"subclass.l1"], #
                           "Cellclass" = data@meta.data$class, #ClusterClass,# 
                           "Patient" = data@meta.data$patient) %>% #orig.ident) %>% # 
      group_by(Disease) %>% dplyr::mutate(Disease = paste0(Disease, " (", n_distinct(Patient), ")")) %>%
      group_by(Celltype) %>% dplyr::mutate(Celltype2 = paste0(Celltype, " (", length(Expression), ")"))
    
    data.sum = data.plot %>%
      group_by(Disease, Celltype2, Cellclass) %>% dplyr::summarize(Percentage = sum(Expression>0)/n(), Expression = mean(Expression))
    
    list("data.plot" = data.plot, "data.sum" = data.sum)
  })
  data.patient = reactive({#eventReactive(c(input$celltype, input$gene, input$pat.grouping, input$cell.grouping),{
    req(input$celltype)
    data.patient = data.celltype()$data.plot %>%
      group_by(Disease, Patient, Celltype) %>% dplyr::summarize(Percentage = sum(Expression>0)/n(), Expression = mean(Expression), .groups = "keep") %>%
      tibble %>% tidyr::complete(Celltype, nesting(Disease, Patient), fill = list("Expression" = 0, "Percentage" = 0)) %>%
      filter(Celltype %in% input$celltype)
    
    data.num = data.celltype()$data.plot %>%
      group_by(Disease, Patient) %>% 
      dplyr::summarize(AllCells = length(Expression),
                       SelCells = length(Expression[Celltype %in% input$celltype])/length(Expression))
    
    list("data.patient" = data.patient, "data.num" = data.num)
  })
  output$plot.celltype = renderPlot({
    data.sum = data.celltype()$data.sum
    dot.plot = ggplot(data.sum, aes(x = Disease, y = Celltype2)) +
      geom_point(aes(size = Percentage, col = Expression)) +
      {if(input$cell.class) facet_grid(Cellclass ~ ., scales = "free_y", space = "free", switch = "y")} +
      scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(data.sum$Percentage)), guide = guide_legend(direction = "vertical")) +
      scale_color_viridis(guide = guide_colorbar(direction = "vertical")) +
      theme_bw() + theme(text = element_text(size = 18), strip.text.y.left = element_text(angle = 0), legend.position = "bottom", 
                         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    data.plot = data.celltype()$data.plot
    box.plot = ggplot(data.plot, aes(x = Expression, y = Celltype2)) +
      geom_boxplot(fill = "lightblue") +
      scale_x_continuous(expand = c(0,0)) + scale_fill_discrete(guide = guide_legend(direction = "vertical")) +
      {if(input$cell.class) facet_grid(Cellclass ~ ., scales = "free", space = "free", switch = "both")} +
      theme_bw() + theme(text = element_text(size = 18), legend.position = "bottom", 
                         axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), 
                         strip.background = element_blank(), strip.text.y = element_blank())
    
    dot.plot + box.plot + plot_layout(widths = c(1,3))
  })
  output$plot.patient = renderPlot({
      data.patient = data.patient()$data.patient
      dot.plot = ggplot(data.patient, aes(x = Patient, y = Celltype)) +
        geom_point(aes(size = Percentage, col = Expression)) +
        facet_grid(. ~ Disease, scales = "free", space = "free", switch = "both") +
        scale_size_continuous(labels = scales::percent, limits = c(0.0000001, max(c(0, data.patient$Percentage), na.rm = TRUE)), guide = guide_legend(direction = "vertical")) +
        scale_color_viridis(guide = guide_colorbar(direction = "vertical")) +
        theme_bw() + theme(text = element_text(size = 18), strip.text.y.left = element_text(angle = 0), legend.position = "bottom", 
                           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      
      data.num = data.patient()$data.num
      bar.plot = ggplot(data.num, aes(x = Patient)) +
        geom_col(aes(y = AllCells), fill = "black", width = 0.4, position = position_nudge(x = -0.2)) +
        geom_col(aes(y = SelCells*max(AllCells)/max(SelCells)), fill = "grey80", width = 0.4, position = position_nudge(x = 0.2)) +
        facet_grid(. ~ Disease, scales = "free", space = "free", switch = "both") +
        scale_y_continuous(name = "All cells", expand = c(0,0), 
                           sec.axis = sec_axis(~./max(data.num$AllCells)*max(data.num$SelCells), name = paste("Selected cell fraction"), labels = scales::percent)) +
        theme_bw() + theme(text = element_text(size = 18), 
                           axis.title.y.right = element_text(color = "grey80"), axis.text.y.right = element_text(color = "grey80"), axis.ticks.y.right = element_line(color = "grey80"),
                           axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), 
                           strip.background = element_blank(), strip.text.x = element_blank())

      bar.plot / dot.plot + plot_layout(heights = c(2,0.75+0.25*length(input$celltype)))
  })
  output$ui.plot.patient = renderUI({
    plotOutput("plot.patient", height = 450+20*length(input$celltype))
  })
  output$download.xlsx = downloadHandler(
    filename = function() { paste0('Excel-data_KPMP-browser_', Sys.Date(), '.xlsx') },
    content = function(file) {
      wb = createWorkbook()
      addDataFrame(as.data.frame(data.celltype()$data.sum), sheet = createSheet(wb, "Data.sum"), row.names=FALSE)
      addDataFrame(as.data.frame(data.patient()$data.patient), sheet = createSheet(wb, "Data.patient"), row.names=FALSE)
      addDataFrame(as.data.frame(data.patient()$data.num), sheet = createSheet(wb, "Data.num"), row.names=FALSE)
      saveWorkbook(wb, file)
  })
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))
