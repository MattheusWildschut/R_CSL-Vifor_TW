## Shiny Image Curator
setwd("C:/Incucyte/Netosis/Netosis_Exp10/")
source("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/SourceFile_TW.R")

## Data Table data.nuc ---------------------------------------------------------------------------------
# crop.size = 50
# data.nuc.raw = fread("CP_Analyses/Netosis-10_Nucleus.csv")
# data.nuc = data.nuc.raw %>%
#   `colnames<-`(str_remove(colnames(.), "Metadata_")) %>%
#   dplyr::select(ImageNumber, ObjectNumber, Experiment, Well, Image, Timepoint, Location_Center_X, Location_Center_Y) %>%
#   mutate(Timepoint = paste0("00d", Timepoint),
#          Filename = paste(Experiment, "Phase", Well, Image, Timepoint, sep = "_"),
#          Filename2 = paste(Well, Image, Timepoint, sep = "_"),
#          X = Location_Center_X-(crop.size/2),
#          Y = Location_Center_Y-(crop.size/2),
#          Border = X < 0 | X > 1264-crop.size | Y < 0 | Y > 936-crop.size,
#          Class = as.character(NA)) %>%
#   filter(!Border & Timepoint != "00d06h00m")
# write.csv(data.nuc, "Classification/Netosis_Exp10_Classification.csv")

## UI ----------------------------------------------------------------------------------------------------------------------------------
ui = fluidPage(
  useShinyjs(),
  extendShinyjs(text = "shinyjs.closewindow = function() { window.close(); }", functions = "closewindow"),
  titlePanel("Shiny Image Curator_TW"),
  sidebarLayout(
    sidebarPanel(
      shinyFilesButton("class", label = "Select class file", title = "Please select class file", multiple = "single",
                       buttonType = "success"),
      verbatimTextOutput("class_path"),
      actionBttn(inputId = "save_close", label = "Save & close", style = "unite", color = "success")
    ),
    mainPanel(
      tags$h3("Draw rectangle & double-click to (reset) zoom | Use class buttons for cell classifications"),
      fluidRow(column(width = 10, imageOutput("im2", height = 745, width = 1000,
                                              dblclick = "im2_dblclick", 
                                              brush = brushOpts(id = "im2_brush", resetOnNew = TRUE)),
                      actionBttn(inputId = "NET", label = "NET", style = "unite", color = "success"),
                      actionBttn(inputId = "SmallGreen", label = "SmallGreen", style = "unite", color = "primary"),
                      actionBttn(inputId = "Original", label = "Original", style = "unite", color = "royal"),
                      actionBttn(inputId = "Adherent", label = "Adherent", style = "unite", color = "danger"),
                      actionBttn(inputId = "FlatCell", label = "FlatCell", style = "unite", color = "warning"),
                      actionBttn(inputId = "NoCell", label = "NoCell", style = "unite", color = "default"),
                      actionBttn(inputId = "previous", label = "Go back to previous cell", style = "unite", color = "danger")),
               column(width = 2, plotOutput("im", height = 745))),
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
    write.csv(data.nuc2(), paste0("Classification/Netosis_Exp10_Classification_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv"),
              row.names = FALSE)
    js$closewindow();
    stopApp()
  })
  roots = c("Input" = "Classification", Home = "C:/Incucyte/Netosis/")
  shinyFileChoose(input, "class", roots=roots, filetypes=c("", "txt", "csv", "xls", "xlsx"))
  class_path = reactive(unlist(parseFilePaths(roots, input$class)[,"datapath"]))
  output$class_path = renderText({
    if(is.integer(input$class)){
      "Please select class file"
    } else {
      paste0("Selected class file:\n", class_path())
    }
  })
  data.nuc2 = reactiveVal(NULL)
  observeEvent(req(!is.integer(input$class)),
    data.nuc2(fread(class_path()))
  )

  buttons = reactive(list(input$NET,input$SmallGreen, input$Original, input$Adherent, input$FlatCell, input$NoCell))
  observeEvent(input$NET, {class("NET")})
  observeEvent(input$SmallGreen, {class("SmallGreen")})
  observeEvent(input$Original, {class("Original")})
  observeEvent(input$Adherent, {class("Adherent")})
  observeEvent(input$FlatCell, {class("FlatCell")})
  observeEvent(input$NoCell, {class("NoCell")})
  
  x = reactiveVal(NULL)
  x_list = reactiveVal(list(NULL))
  brush = reactiveVal(NULL)
  
  observeEvent(buttons(), {
    data.nuc2(data.nuc2() %>%
                mutate(Class = ifelse(ObjectNumber == obj.n() & ImageNumber == img.n(), class(), Class)))
    x_list(append(x_list(), x()))
    brush(NULL)
  }, ignoreInit = TRUE)
  observeEvent(data.nuc2(), {
    x(sample(which(is.na(data.nuc2()$Class)),1))
  })
  observeEvent(input$previous, {
    x(x_list()[[length(x_list())]])
    x_list(x_list()[1:(length(x_list())-1)])
  })

  img.n = reactiveVal(NULL)
  obj.n = reactiveVal(NULL)
  class = reactiveVal(NULL)
  observeEvent(x(), {
    img.n(data.nuc2()$ImageNumber[x()])
    obj.n(data.nuc2()$ObjectNumber[x()])
  })
  
  observeEvent(input$im2_dblclick, {
    brush(input$im2_brush)
  })
  output$im2 = renderImage({
    req(x())
    img = image_read(paste0("MultiChannel_Images/Netosis_Exp10_Merged_", data.nuc$Filename2[x()], ".jpg"))
    img2 = image_draw(img)
    coord = data.nuc[x(), c("X", "Y")]
    rect(coord$X, coord$Y, coord$X+50, coord$Y+50, border = "white", lty = "dashed", lwd = 3)
    img2 = image_resize(img2, "1000x740")
    if(!is.null(brush())){
      img2 = image_crop(img2, geometry_area(width = brush()$xmax-brush()$xmin, height = brush()$ymax-brush()$ymin, 
                                            x_off = brush()$xmin, y_off = brush()$ymin))
    }
    outfile <- tempfile(fileext='.png')
    image_write(image_resize(img2, "1000x740"), outfile)
    dev.off()
    list(src = outfile, width = "auto", height = "auto")
  }, deleteFile = TRUE)
  output$im = renderPlot({
    req(x())
    im = image_read(paste0("Cell_Crops/Netosis-10_Phase_", data.nuc$Filename2[x()], "_", obj.n(), ".tif"))[1]
    im2 = image_read(paste0("Cell_Crops/Netosis-10_Red_", data.nuc$Filename2[x()], "_", obj.n(), ".tif"))[1]
    im3 = image_read(paste0("Cell_Crops/Netosis-10_Green_", data.nuc$Filename2[x()], "_", obj.n(), ".tif"))[1]
    im = image_level(im, black_point = input$Phase[1], white_point = input$Phase[2], mid_point = input$Phase2)
    im2 = image_level(im2, black_point = input$Red[1], white_point = input$Red[2], mid_point = input$Red2)
    im3 = image_level(im3, black_point = input$Green[1], white_point = input$Green[2], mid_point = input$Green2)
    image_ggplot(im) / image_ggplot(im2) / image_ggplot(im3)
  })
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))
