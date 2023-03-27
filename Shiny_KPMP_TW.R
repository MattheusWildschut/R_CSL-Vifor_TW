setwd("C:\\Users\\mwildschut\\OneDrive - Vifor Pharma AG\\Documents\\R\\Projects\\R_CSL-Vifor_TW")
source("SourceFile_TW.R")

# install.packages("Seurat")
# remotes::install_github("mojaveazure/seurat-disk")
# install.packages("Matrix")

library(Seurat)
library(SeuratDisk)

# data.scRNA = LoadH5Seurat("Input\\GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat", assays = list(RNA = "data")) #%>%
  # subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
# data.snRNA = LoadH5Seurat("Input\\GSE183277_Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat", assays = list(RNA = "data"))# %>%
  # subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
annot.scRNA = fread("Input//GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_Field_Descriptions.txt.gz") 
annot.snRNA = fread("Input//GSE183277_Kidney_Healthy-Injury_Cell_Atlas_snCv3_Metadata_Field_Descriptions.txt.gz")
# fields.pat = data.frame("Field" = c(paste0("condition", c(".long", ".l1", ".l2", ".l3")), "patient", "sex", "race")) %>%
#   mutate(Description = annot.scRNA$Description[match(Field, annot.scRNA$Field)])
# fields.cell = data.frame("Field" = c("class", paste0("subclass", c(".l1", ".l2", ".l3", ".full")), "state.l1", "state.l2")) %>%
#   mutate(Description = annot.scRNA$Description[match(Field, annot.scRNA$Field)])

## Data from https://atlas.kpmp.org/repository/?facetTab=files&files_size=60&files_sort=%5B%7B%22field%22%3A%22
## file_name%22%2C%22order%22%3A%22asc%22%7D%5D&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in
## %22%2C%22content%22%3A%7B%22field%22%3A%22dois%22%2C%22value%22%3A%5B%2210.48698%2F92nk-e805%22%5D%7D%7D%5D%7D
data.scRNA = LoadH5Seurat("Input/521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat", assays = list(RNA = c("data", "counts")))
# data.scRNA = data.scRNA %>% subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 50)
data.snRNA = LoadH5Seurat("Input/c798e11b-bbde-45dd-bd91-487f27c93f8f_WashU-UCSD_HuBMAP_KPMP-Biopsy_10X-R_12032021.h5Seurat", assays = list(RNA = c("data", "counts")))
# data.snRNA = data.snRNA %>% subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 50)
fields.pat = data.frame("Field" = c("sampletype", "diseasetype", "age", "gender", "state", "tissuetype", "celltype")) %>%
  mutate(Description = annot.scRNA$Description[match(Field, annot.scRNA$Field)])
fields.cell = data.frame("Field" = c("ClusterClass", paste0("subclass", c(".l1", ".l2")), "celltype")) %>%
  mutate(Description = annot.scRNA$Description[match(Field, annot.scRNA$Field)])

## UI ----------------------------------------------------------------------------------------------------------------------------------
ui = {fluidPage(
  useShinyjs(),
  titlePanel("CSL-Vifor_KPMP-browser_TW"),
  sidebarLayout(
    sidebarPanel(
                 radioGroupButtons(inputId = "dataset",
                                   label = "Choose scRNA-seq or snRNA-seq dataset",
                                   choices = c("scRNA", "snRNA"),
                                   justified = TRUE,
                                   checkIcon = list(yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                                                    no = tags$i(class = "fa fa-circle-o", style = "color: steelblue"))),
                 selectizeInput(inputId = "gene",
                                label = "Select gene of interest",
                                choices = NULL),
                 pickerInput(inputId = "pat.grouping",
                             label = "Select patient grouping", 
                             choices = fields.pat$Field,
                             # choicesOpt = list(subtext = fields.pat$Description),
                             selected = "diseasetype"), #"condition.l1"
                 pickerInput(inputId = "cell.grouping",
                             label = "Select cell grouping", 
                             choices = fields.cell$Field,
                             # choicesOpt = list(subtext = fields.cell$Description),
                             selected = "subclass.l1"),
                 materialSwitch(inputId = "cell.class",
                                label = "Split cell types by class",
                                value = TRUE,
                                status = "primary"),
                 uiOutput(outputId = "ui.celltype"),
                 pickerInput(inputId = "umap.celltypes",
                             label = "Choose celltypes included in UMAP",
                             choices = sort(as.character(unique(data.scRNA@meta.data$subclass.l1))),
                             selected = sort(as.character(unique(data.scRNA@meta.data$subclass.l1))),
                             multiple = TRUE,
                             options = list(`actions-box` = TRUE)),
                 downloadBttn(outputId = "download.xlsx",
                              label = list(icon("file-excel"), "Download data"),
                              color = "success", size = "sm"),
                 downloadBttn(outputId = "download.celltype",
                              label = list(icon("file-pdf"), "Download Plot 1"),
                              color = "danger", size = "sm"),
                 downloadBttn(outputId = "download.patient",
                              label = list(icon("file-pdf"), "Download Plot 2"),
                              color = "danger", size = "sm")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("BubblePlots",
                 plotOutput("plot.celltype", height = 1000),
                 uiOutput("ui.plot.patient")),
        tabPanel("UMAP", 
                 div(
                   div(style="display: inline-block; width: 700px; vertical-align:top", 
                       uiOutput("ui.plot.umap1")),
                   div(style="display: inline-block; width: 525px; vertical-align:top", 
                       plotOutput("plot.umap2", height = 500, width = 525))
                 )
        )
      )
    )
  )
)}

## Server ----------------------------------------------------------------------------------------------------------------------------------
server = function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })
  data = reactive({
    # get(paste0("data.", input$dataset))
    data.scRNA
  })
  # Gene expression Server ---------
  updateSelectizeInput(session, "gene", choices = sort(rownames(data.scRNA@assays$RNA@data)), 
                       selected = "P2RY14", options = list(maxOptions = 100000), server = TRUE)
  output$ui.celltype = renderUI({
    celltypes = unique(sort(as.character(unlist(data()@meta.data[,input$cell.grouping]))))
    pickerInput(inputId = "celltype",
                label = "Select cell types of interest for single patient assesment",
                choices = celltypes,
                multiple = TRUE,
                selected = NULL,
                options = list(`actions-box` = TRUE))
  })
  data.celltype = reactive({
    req(input$gene, input$pat.grouping, input$cell.grouping)
    data.plot = data.frame("Expression" = data()@assays$RNA@data[input$gene,], #["SLC9B2",], #
                           "Disease" = data()@meta.data[,input$pat.grouping], #[,"condition.l1"], #
                           "Celltype" = data()@meta.data[,input$cell.grouping], #[,"subclass.l1"], #
                           "Cellclass" = data()@meta.data$ClusterClass, #ClusterClass,# 
                           "Patient" = data()@meta.data$SpecimenID) %>% #orig.ident) %>% # 
      group_by(Disease) %>% dplyr::mutate(Disease = paste0(Disease, " (", n_distinct(Patient), ")")) %>%
      group_by(Celltype) %>% dplyr::mutate(Celltype2 = paste0(Celltype, " (", length(Expression), ")"))
    
    data.sum = data.plot %>%
      group_by(Disease, Celltype2, Cellclass) %>% dplyr::summarize(Percentage = sum(Expression>0)/n(), Expression = mean(Expression))
    
    list("data.plot" = data.plot, "data.sum" = data.sum)
  })
  data.patient = reactive({
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
  # data.umap = reactive({
  #   req(input$gene)
  #   data.plot = data.frame("UMAP_1" = round(data.scRNA@meta.data$umap1, 2),
  #                          "UMAP_2" = round(data.scRNA@meta.data$umap2, 2),
  #                          "CellType" = data.scRNA@meta.data$subclass.l1,
  #                          "CellType2" = data.scRNA@meta.data$subclass.l2,
  #                          "Expression" = data.scRNA@assays$RNA@data[input$gene,])
  # })
  plot.celltype = reactive({
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
  plot.patient = reactive({
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
  # output$plot.umap1 = renderPlot({
  output$plot.umap1 = renderPlot({
    # ggplot(data.umap(), aes(x = UMAP_1, UMAP_2, col = CellType)) +
    #   geom_point(size = 0.1) + theme_bw()
    DimPlot(data.scRNA, reduction = "umap", group.by = input$cell.grouping)
  })
  output$ui.plot.umap1 = renderUI({
    plotOutput("plot.umap1", height = 500, 
               width = switch(input$cell.grouping,
                              "ClusterClass" = 600,
                              "subclass.l1" = 575,
                              700))
  })
  output$plot.umap2 = renderPlot({
    req(input$umap.celltypes, input$gene)
    FeaturePlot(data.scRNA, reduction = "umap", features = input$gene, cols = c("grey95", viridis(3)),
                cells = WhichCells(data.scRNA, expression = subclass.l1 %in% input$umap.celltypes))
    # data.plot = filter(data.umap(), CellType %in% input$umap.celltypes)
    # ggplot(data.plot, aes(x = UMAP_1, UMAP_2, col = Expression,
    #                         text = paste0("Class: ", CellType, "\nSubtype: ", CellType2))) +
    #   geom_point(size = 0.1) +
    #   scale_color_viridis(limits = c(0.00001, max(c(0.00001, data.plot$Expression))), na.value = "grey90") + 
    #   theme_bw()
  })
  # output$ui.plot.umap2 = renderUI({
  #   plotOutput("plot.umap2", height = 500)
  # })
  output$plot.celltype = renderPlot({
    plot.celltype()
  })
  output$plot.patient = renderPlot({
    plot.patient()
  })
  output$ui.plot.patient = renderUI({
    plotOutput("plot.patient", height = 450+20*length(input$celltype))
  })
  output$download.xlsx = downloadHandler(
    filename = function() { paste0('Excel-data_KPMP-browser_', Sys.Date(), '.xlsx') },
    content = function(file) {
      wb = createWorkbook()
      addDataFrame(as.data.frame(data.celltype()$data.sum %>% mutate(Percentage = percent(Percentage))),
                   sheet = createSheet(wb, "Data.sum"), row.names = FALSE)
      if(!is.null(input$celltype)){
        addDataFrame(as.data.frame(data.patient()$data.patient %>% mutate(Percentage = percent(Percentage))), 
                     sheet = createSheet(wb, "Data.patient"), row.names = FALSE)
        addDataFrame(as.data.frame(data.patient()$data.num %>% mutate(SelCells = percent(SelCells))),
                                   sheet = createSheet(wb, "Data.num"), row.names = FALSE)
      }
      saveWorkbook(wb, file)
  })
  output$download.celltype = downloadHandler(
    filename = function() { paste0("Celltype-plot_", input$gene, "_", input$pat.grouping, "-", input$cell.grouping, "_", Sys.Date(), '.pdf') },
    content = function(file) {
      ggsave(file, plot.celltype(), height = 16, width = 9)
    })
  output$download.patient = downloadHandler(
    filename = function() { paste0("Patient-plot_", input$gene, "_", input$pat.grouping, "-", input$cell.grouping, "_", Sys.Date(), '.pdf') },
    content = function(file) {
      if(!is.null(input$celltype)){
        ggsave(file, plot.patient(), height = 8, width = 10)  
      } else {
        showNotification("Select cell type of interest first")
      }
    }
    )
  observe({
    if(is.null(input$celltype)) shinyjs::hide("download.patient_bttn") else shinyjs::show("download.patient_bttn")
  })
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))


data.scRNA.raw = LoadH5Seurat("Input\\521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat", assays = list(RNA = c("data", "counts")))
data.scRNA = data.scRNA.raw %>% subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)

d = AverageExpression(data.scRNA, group.by = c("subclass.l2", "diseasetype"))[[1]]

data.plot = data.frame("Expression" = data.scRNA@assays$RNA@data["P2RY14",], #["SLC9B2",], #
                       "Disease" = data.scRNA@meta.data$diseasetype, #[,"condition.l1"], #
                       "Celltype" = data.scRNA@meta.data$subclass.l2, #[,"subclass.l1"], #
                       "Cellclass" = as.character(data.scRNA@meta.data$subclass.l1), #ClusterClass,# 
                       "Patient" = data.scRNA@meta.data$SpecimenID) %>% #orig.ident) %>% # 
  group_by(Disease) %>% dplyr::mutate(Disease2 = paste0(Disease, " (", n_distinct(Patient), ")")) %>%
  group_by(Celltype) %>% dplyr::mutate(Celltype2 = paste0(Celltype, " (", length(Expression), ")"))# %>%
  # filter(Celltype2 %in% c("cDC (1708)", "pDC (207)"))
  # filter(Celltype2 == "pDC (207)" & Disease == "LivingDonor (20)") %>%
  # mutate(LOG = 2^Expression)

data.sum = data.plot %>%
  group_by(Disease, Disease2, Celltype, Celltype2, Cellclass) %>% 
  dplyr::summarize(Count = n(), Percentage = sum(Expression>0)/n())

data.avg = data.frame("Expression" = d["P2RY14",],
                       "Disease" = str_split_fixed(colnames(d), "_", 2)[,2],
                       "Celltype" = str_split_fixed(colnames(d), "_", 2)[,1])

data.sum2 = left_join(data.sum, data.avg, by = c("Disease", "Celltype")) %>%
  # rbind(data.sum %>% mutate(Percentage = Percentage*n))
  mutate(Count2 = Count/as.numeric(str_extract(Celltype2, "(?<=\\()[:digit:]+")))


ggplot(data.sum2, aes(x = Disease2, y = Celltype2)) +
  geom_tile(aes(fill = Count2), width = 0.4, col = "white") +
  scale_fill_gradient(name = "Percent cells", low = "white", high = "grey40", guide = guide_colorbar(direction = "vertical"), labels = scales::percent, limits = c(0.0000001, max(data.sum2$Count2))) +
  geom_point(aes(col = Expression, size = Percentage)) + #size = Percentage, 
  scale_y_discrete(limits = rev) +
  facet_grid(Cellclass ~ ., scales = "free_y", space = "free", switch = "y") +
  scale_size_continuous(name = "Percent expressing", guide = guide_legend(direction = "vertical"), 
                        labels = scales::percent, limits = c(0.0000001, max(data.sum2$Percentage))) +
  scale_color_viridis(guide = guide_colorbar(direction = "vertical")) +
  xlab(NULL) + ylab("Cell type") +
  # scale_color_continuous(guide = guide_colorbar(direction = "vertical", reverse = TRUE), trans = "reverse") +
  # scale_fill_distiller(trans = "reverse", guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  # scale_color_gradient2(midpoint = mean(data.sum$Count)) +#palette = "YlOrBr", trans = "reverse", guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  theme_bw() + theme(text = element_text(size = 18), strip.text.y.left = element_text(angle = 0), legend.position = "bottom", 
                     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

DimPlot(data.scRNA, reduction = "umap", label = TRUE, group.by = "subclass.l1")
a = data.scRNA@meta.data
data.plot = data.frame("UMAP_1" = round(data.scRNA@meta.data$umap1, 2),
                       "UMAP_2" = round(data.scRNA@meta.data$umap2, 2),
                       "CellType" = data.scRNA@meta.data$subclass.l1,
                       "CellType2" = data.scRNA@meta.data$subclass.l2,
                       "Disease" = data.scRNA@meta.data$diseasetype,
                       "Expression" = data.scRNA@assays$RNA@data["P2RY14",])
ggplotly(
  ggplot(data.plot, aes(x = UMAP_1, UMAP_2, col = Expression, shape = Disease, fill = CellType,
                        text = paste0("Class: ", CellType, "\nSubtype: ", CellType2))) +
    geom_point(size = 0.1) + 
    scale_color_viridis(limits = c(0.00001, max(data.plot$Expression)), na.value = "grey90") + theme_bw()
)

ggplot(data.plot, aes(x = UMAP_1, UMAP_2, col = CellType)) +
  geom_point(size = 0.1) + theme_bw()


DimPlot(data.scRNA, reduction = "umap", group.by = "subclass.l1")
FeaturePlot(data.scRNA, reduction = "umap", features = "P2RY14", cols = c("grey95", viridis(3)),
            cells = WhichCells(data.scRNA, expression = subclass.l1=="Immune"))

as.formula("a")
