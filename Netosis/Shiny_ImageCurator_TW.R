## Shiny Image Curator
setwd("C:/Incucyte/Netosis/")
source("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/SourceFile_TW.R")
library(shinyFiles)
library(magick)
library(shinyjs)

## Data Table data.nuc ---------------------------------------------------------------------------------
exp.folder = "C:/Incucyte/Netosis/Netosis_Exp10/CP_Analyses/"
crop.size = 50
data.nuc.raw = fread(paste0(exp.folder, "Netosis-10_Nucleus.csv"))
data.nuc = data.nuc.raw %>%
  `colnames<-`(str_remove(colnames(.), "Metadata_")) %>%
  dplyr::select(ImageNumber, ObjectNumber, Experiment, Well, Image, Timepoint, Location_Center_X, Location_Center_Y) %>%
  mutate(Timepoint = paste0("00d", Timepoint),
         Filename = paste(Experiment, "Phase", Well, Image, Timepoint, sep = "_"),
         Filename2 = paste(Well, Image, Timepoint, sep = "_"),
         X = Location_Center_X-(crop.size/2),
         Y = Location_Center_Y-(crop.size/2),
         Border = X < 0 | X > 1264-crop.size | Y < 0 | Y > 936-crop.size,
         Class = as.character(NA)) %>%
  filter(!Border)

## UI ----------------------------------------------------------------------------------------------------------------------------------
ui = fluidPage(
  useShinyjs(),
  extendShinyjs(text = "shinyjs.closewindow = function() { window.close(); }", functions = "closewindow"),
  titlePanel("Shiny Image Curator_TW"),
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("data_folder", label = "Select folder", title = "Please select folder", buttonType = "success"),
      verbatimTextOutput("folder_path"),
      actionBttn(inputId = "save_close", label = "Save & close", style = "unite", color = "success")
    ),
    mainPanel(
      tags$h3("Use buttons to classify cells"),
      fluidRow(column(width = 6, plotOutput("im"),
                      actionBttn(inputId = "NET", label = "NET", style = "unite", color = "success"),
                      actionBttn(inputId = "SmallGreen", label = "SmallGreen", style = "unite", color = "primary"),
                      actionBttn(inputId = "Original", label = "Original", style = "unite", color = "royal"),
                      actionBttn(inputId = "FlatCell", label = "FlatCell", style = "unite", color = "danger"),
                      actionBttn(inputId = "NoCell", label = "NoCell", style = "unite", color = "warning")),
               column(width = 6, imageOutput("im2"))),
      hr(style = "border-color: black; color: black;"),
      tags$h3("Adjust display of cell crop images"),
      fluidRow(column(width = 4,
                      sliderInput(inputId = "Phase", label = "Phase min/max threshold", min = 0, max = 100, value = c(0,100), step = 1),
                      sliderInput(inputId = "Phase2", label = "Phase midpoint", min = 0, max = 2, value = 1, step = 0.01)),
               column(width = 4,
                      sliderInput(inputId = "Red", label = "Red min/max threshold", min = 0, max = 2, value = c(0,0.75), step = 0.01),
                      sliderInput(inputId = "Red2", label = "Red midpoint", min = 0, max = 2, value = 0.5, step = 0.01)),
               column(width = 4,
                      sliderInput(inputId = "Green", label = "Green min/max threshold", min = 0, max = 25, value = c(0,2), step = 0.01),
                      sliderInput(inputId = "Green2", label = "Green midpoint", min = 0, max = 2, value = 0.5, step = 0.01)))
      )
  )
)

## Server ----------------------------------------------------------------------------------------------------------------------------------
server = function(input, output, session) {
  session$onSessionEnded(function() {stopApp()})
  observeEvent(input$save_close, {
    write.csv(data.nuc2(), "C:/Incucyte/Netosis/Netosis_Exp10/Netosis_Exp10_Classification.csv")
    js$closewindow();
    stopApp()
  })
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
  data.nuc2 = reactiveVal(data.nuc)
  x = reactiveVal(sample(which(is.na(data.nuc$Class)),1))
  img.n = reactive(data.nuc$ImageNumber[x()])
  obj.n = reactive(data.nuc$ObjectNumber[x()])
  class = reactiveVal(NA)
  
  buttons = reactive(list(input$NET,input$SmallGreen, input$Original, input$FlatCell, input$NoCell))
  observeEvent(input$NET, {class("NET")})
  observeEvent(input$SmallGreen, {class("SmallGreen")})
  observeEvent(input$Original, {class("Original")})
  observeEvent(input$FlatCell, {class("FlatCell")})
  observeEvent(input$NoCell, {class("NoCell")})
  
  observeEvent(buttons(), {
    data.nuc = data.nuc2()
    data.nuc2(data.nuc %>%
                mutate(Class = ifelse(ObjectNumber == obj.n() & ImageNumber == img.n(), class(), Class)))
    x(sample(which(is.na(data.nuc2()$Class)),1))
  })
  output$img.n = renderText(img.n())
  output$obj.n = renderText(obj.n())
  output$tab = renderDataTable({
    data.nuc2()
  })
  output$value = renderText(x())
  output$im = renderPlot({
    folder = "C:/Incucyte/Netosis/Netosis_Exp10/Cell_Crops"
    number = x()
    im = image_read(paste0(folder, "/Netosis-10_Phase_", data.nuc$Filename2[x()], "_", obj.n(), ".tif"))[1]
    im2 = image_read(paste0(folder, "/Netosis-10_Red_", data.nuc$Filename2[x()], "_", obj.n(), ".tif"))[1]
    im3 = image_read(paste0(folder, "/Netosis-10_Green_", data.nuc$Filename2[x()], "_", obj.n(), ".tif"))[1]
    im = image_level(im, black_point = input$Phase[1], white_point = input$Phase[2], mid_point = input$Phase2)
    im2 = image_level(im2, black_point = input$Red[1], white_point = input$Red[2], mid_point = input$Red2)
    im3 = image_level(im3, black_point = input$Green[1], white_point = input$Green[2], mid_point = input$Green2)
    image_ggplot(im) + image_ggplot(im2) + image_ggplot(im3)
  })
  output$im2 = renderImage({
    img = image_read(paste0("C:/Incucyte/Netosis/Netosis_Exp10/MultiChannel_Images/Netosis_Exp10_Merged_", data.nuc$Filename2[x()], ".jpg"))
    img2 = image_draw(img)
    coord = data.nuc[x(), c("X", "Y")]
    rect(coord$X, coord$Y, coord$X+50, coord$Y+50, border = "white", lty = "dashed", lwd = 3)
    outfile <- tempfile(fileext='.png')
    image_write(img2, outfile)
    dev.off()
    list(src = outfile, width = 540,  height = 400)
  }, deleteFile = TRUE)
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))
