## Shiny Image Curator
setwd("C:/Incucyte/Netosis/")
source("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/SourceFile_TW.R")
library(shinyFiles)
library(magick)

## UI ----------------------------------------------------------------------------------------------------------------------------------
ui = fluidPage(
  titlePanel("Shiny Image Curator_TW"),
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("data_folder", label = "Select folder", title = "Please select folder", buttonType = "success"),
      verbatimTextOutput("folder_path")
    ),
    mainPanel(
      plotOutput("im"),
      # imageOutput("im2"),
      textOutput("value"),
      column(width = 4,
             sliderInput(inputId = "Phase", label = "Phase min/max threshold", min = 0, max = 100, value = c(0,100), step = 1),
             sliderInput(inputId = "Phase2", label = "Phase midpoint", min = 0, max = 2, value = 1, step = 0.01)),
      column(width = 4,
             sliderInput(inputId = "Red", label = "Red min/max threshold", min = 0, max = 2, value = c(0,0.75), step = 0.01),
             sliderInput(inputId = "Red2", label = "Red midpoint", min = 0, max = 2, value = 0.5, step = 0.01)),
      column(width = 4,
             sliderInput(inputId = "Green", label = "Green min/max threshold", min = 0, max = 25, value = c(0,2), step = 0.01),
             sliderInput(inputId = "Green2", label = "Green midpoint", min = 0, max = 2, value = 0.5, step = 0.01)),
      actionBttn(inputId = "Next", label = "Next image", style = "unite", color = "succes"),
      actionBttn(inputId = "Previous", label = "Previous image", style = "unite", color = "danger")
    #   tabsetPanel(
    #     tabPanel("Data table", dataTableOutput("data"),
    #              uiOutput("ui.download.xlsx")),
    #     tabPanel("Plot | Influence medium", uiOutput("medium.ui"), uiOutput("medium.ui2")),
    #     tabPanel("Plot | qPCR data", uiOutput("plot.ui"), uiOutput("plot.ui2"))
    #   )
    )
  )
)

## Server ----------------------------------------------------------------------------------------------------------------------------------
server = function(input, output, session) {
  session$onSessionEnded(function() {stopApp()})
  roots = c("Input" = "C:/Incucyte/Netosis/Netosis_Exp10/", Home = "C:/Incucyte/Netosis/")
  shinyDirChoose(input, "data_folder", roots=roots, filetypes=c("", "txt", "csv", "xls", "xlsx"))
  # folder <- reactive(unlist(parseDirPath(roots, input$data_folder)))
  
  output$folder_path = renderText({
    if(is.integer(input$data_folder)){
      "Please select data folder"
    } else {
      paste0("Selected data folder:\n", folder())
    }
  })
  x = reactiveVal(1)
  observeEvent(input$Next, {
    x(x()+1)
  })
  observeEvent(input$Previous, {
    x(x()-1)
  })
  output$value = renderText(x())
  output$im = shiny::renderPlot({
    folder = "C:/Incucyte/Netosis/Netosis_Exp10/Cell_Crops"
    number = x()
    im = image_read(paste0(folder, "/Netosis-10_Phase_C10_1_00d08h00m_", x(), ".tif"))[1] #filename, ".tif"))[1]
    im2 = image_read(paste0(folder, "/Netosis-10_Red_C10_1_00d08h00m_", x(), ".tif"))[1]
    im3 = image_read(paste0(folder, "/Netosis-10_Green_C10_1_00d08h00m_", x(), ".tif"))[1]
    im = image_level(im, black_point = input$Phase[1], white_point = input$Phase[2], mid_point = input$Phase2)
    im2 = image_level(im2, black_point = input$Red[1], white_point = input$Red[2], mid_point = input$Red2)
    im3 = image_level(im3, black_point = input$Green[1], white_point = input$Green[2], mid_point = input$Green2)
    image_ggplot(im) + image_ggplot(im2) + image_ggplot(im3)
  })
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))

im2 = image_read(paste0(folder, "/Netosis-10_Red_C10_1_00d08h00m_", 1, ".tif"))[1]
image_normalize(im)
im2 = image_level(im2, black_point = 0, white_point = .75, mid_point = .5, channel = "Red")
image_ggplot(im2, interpolate = TRUE)
image_level(im2, black_point = 0, white_point = 2, mid_point = .5, channel = "Red")

image_colorize(im2, 20, "red")
image_fill(im2, "red")

a = image_ggplot(im2)
image_ggplot(im2) + scale_fill_viridis()
a$scales$scales

im2j = im2 %>% image_write(tempfile(fileext='.jpg'), format = 'jpg')
image_read(im2j, channel = "red")

image_draw(im2, pointsize = 2, res = .5)
a = as.data.frame(as.integer(im2[[1]]))
col_fun = colorRamp2(c(min(a), max(a)), c("black", "red"))
red = as_ggplot(grid.grabExpr(draw(Heatmap(as.matrix(as.data.frame(as.integer(im2[[1]]))), col = col_fun, 
                                           cluster_rows = FALSE, cluster_columns = FALSE,
                                           show_column_names = FALSE, show_heatmap_legend = FALSE))))

as_ggplot(grid.grabExpr(draw(Heatmap(matrix(1:10, ncol =2)))))
dev.off()
dev.new()

library(EBImage)
img = readImage(paste0(folder, "/Netosis-10_Red_C10_1_00d08h00m_", 1, ".tif"))
plot_ly(x = 1:50, y = 1:50, z = img@.Data, type="heatmap", colors = colorRamp(c("black", "red")))
img = readImage('https://upload.wikimedia.org/wikipedia/commons/thumb/0/00/Crab_Nebula.jpg/240px-Crab_Nebula.jpg')
volcano
