## Shiny Image Curator
source("C:/Users/mwildschut/OneDrive - CSL Behring/Documents/R/Projects/R_CSL-Vifor_TW/SourceFile_TW.R")

## UI ----------------------------------------------------------------------------------------------------------------------------------
ui = fluidPage(
  useShinyjs(),
  extendShinyjs(text = "shinyjs.closewindow = function() { window.close(); }", functions = "closewindow"),
  titlePanel("Shiny Image Curator_TW"),
  tabsetPanel(
    tabPanel("Curation",
             sidebarLayout(
               sidebarPanel(
                 radioGroupButtons(inputId = "exp",
                                   label = "Choose experiment",
                                   choices = str_remove(list.files("C:/Incucyte/Netosis", pattern = "Netosis_Exp"), "Netosis_"),
                                   selected = character(0),
                                   justified = TRUE,
                                   checkIcon = list(yes = icon("circle", list("fa-solid"), style = "color: steelblue"),
                                                    no = icon("circle", style = "color: steelblue"))),
                 shinyFilesButton("class_file", label = "Select class file", title = "Please select class file", multiple = "single",
                                  buttonType = "success"),
                 verbatimTextOutput("class_path"),
                 hr(style = "border-color: black; color: black;"),
                 tableOutput("n_class1"),
                 hr(style = "border-color: black; color: black;"),
                 # verbatimTextOutput("x"),
                 # verbatimTextOutput("img_n"),
                 # verbatimTextOutput("obj_n"),
                 verbatimTextOutput("img_name"),
                 searchInput(
                   inputId = "crop_n",
                   label = "Select specific cell crop to correct", 
                   placeholder = "Enter crop number (see 'Classified')",
                   btnSearch = icon("search"), btnReset = icon("remove")
                 ),
                 hr(style = "border-color: black; color: black;"),
                 actionBttn(inputId = "save_close", style = "unite", color = "success", 
                            label = list(icon("floppy-disk", list("fa-solid", "fa-lg")), "Save & close", icon("circle-xmark", list("fa-solid", "fa-lg"))))
                 ),
               mainPanel(
                 tags$h3("Draw rectangle & double-click to (reset) zoom | Use class buttons for cell classifications"),
                 fluidRow(column(width = 10, imageOutput("im2", height = 745, width = 1000,
                                                         dblclick = "im2_dblclick", 
                                                         brush = brushOpts(id = "im2_brush", resetOnNew = TRUE)),
                                 actionBttn(inputId = "NET", label = "NET", style = "unite", color = "success"),
                                 # actionBttn(inputId = "SmallGreen", label = "SmallGreen", style = "unite", color = "primary"),
                                 actionBttn(inputId = "Round", label = "Round", style = "unite", color = "primary"),
                                 actionBttn(inputId = "Adherent", label = "Adherent", style = "unite", color = "danger"),
                                 actionBttn(inputId = "Flat", label = "Flat", style = "unite", color = "warning"),
                                 actionBttn(inputId = "NoCell", label = list(icon("xmark"), "NoCell"), style = "unite", color = "royal"),
                                 actionBttn(inputId = "previous", label = list(icon("arrow-rotate-left"), "Go back to previous cell"), style = "unite", color = "default")))#,
                          # column(width = 2, plotOutput("im", height = 745))),
                 # hr(style = "border-color: black; color: black;"),
                 # tags$h3("Adjust display of cell crop images"),
                 # fluidRow(column(width = 4,
                 #                 sliderInput(inputId = "Phase", label = "Phase min/max threshold", min = 0, max = 100, value = c(0,100), step = 1),
                 #                 sliderInput(inputId = "Phase2", label = "Phase midpoint", min = 0, max = 2, value = 1, step = 0.01)),
                 #          column(width = 4,
                 #                 sliderInput(inputId = "Red", label = "Red min/max threshold", min = 0, max = 2, value = c(0,0.75), step = 0.01),
                 #                 sliderInput(inputId = "Red2", label = "Red midpoint", min = 0, max = 2, value = 0.5, step = 0.01)),
                 #          column(width = 4,
                 #                 sliderInput(inputId = "Green", label = "Green min/max threshold", min = 0, max = 25, value = c(0,5), step = 0.01),
                 #                 sliderInput(inputId = "Green2", label = "Green midpoint", min = 0, max = 2, value = 1, step = 0.01)))
               )
             )
    ),
    tabPanel("Classified",
             sidebarLayout(
               sidebarPanel(
                 radioGroupButtons(inputId = "class",
                                   label = "Choose class to display classified cell crops",
                                   choices = c("NET", "Round", "Adherent", "Flat", "NoCell"),
                                   justified = TRUE,
                                   checkIcon = list(yes = icon("circle", list("fa-solid"), style = "color: steelblue"),
                                                    no = icon("circle", style = "color: steelblue"))),
                 tableOutput("n_class2")
               ),
               mainPanel(
                 plotOutput("img_class", inline = FALSE, fill = FALSE, width = 1200, height = 720)
               )
             )
             
    )
  )
)

## Server ----------------------------------------------------------------------------------------------------------------------------------
server = function(input, output, session) {
  session$onSessionEnded(function() {stopApp()})
  observe({
    req(input$exp)
    setwd(paste0("C:/Incucyte/Netosis/Netosis_", input$exp))
  })
  roots = c("Input" = "Classification", "CellType" = "CNN_CellType/Classification", Home = "C:/Incucyte/Netosis/")
  shinyFileChoose(input, "class_file", roots=roots, filetypes=c("", "txt", "csv", "xls", "xlsx"))
  class_path = reactive(unlist(parseFilePaths(roots, input$class_file)[,"datapath"]))
  output$class_path = renderText({
    if(is.integer(input$class_file)){
      "Please select class file"
    } else {
      paste0("Selected class file:\n", class_path())
    }
  })
  observeEvent(input$save_close, {
    write.csv(data.nuc2(), str_replace(class_path(), "_2023.*.csv|.csv", paste0("_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv")),
              row.names = FALSE)
    js$closewindow();
    stopApp()
  })
  data.nuc2 = reactiveVal(NULL)
  observeEvent(req(!is.integer(input$class_file)),
    data.nuc2(fread(class_path()) %>%
      dplyr::mutate(Class = factor(Class, levels = c("NET", "Round", "Adherent", "Flat", "NoCell"))))
  )

  buttons = reactive(list(input$NET, input$Round, input$Adherent, input$Flat, input$NoCell))
  observeEvent(input$NET, {class("NET")})
  observeEvent(input$Round, {class("Round")})
  observeEvent(input$Adherent, {class("Adherent")})
  observeEvent(input$Flat, {class("Flat")})
  observeEvent(input$NoCell, {class("NoCell")})
  
  x = reactiveVal(NULL)
  x_list = reactiveVal(list(NULL))
  brush = reactiveVal(NULL)
  
  observeEvent(buttons(), {
    data.nuc2(data.nuc2() %>%
                mutate(Class = ifelse(Filename2 == Filename2[x()] & ObjectNumber == ObjectNumber[x()], class(), as.character(Class))) %>%
                mutate(Class = factor(Class, levels = c("NET", "Round", "Adherent", "Flat", "NoCell")))
              )
    x_list(append(x_list(), x()))
    brush(NULL)
    updateSearchInput(session = session, inputId = "crop_n", value = "")
  }, ignoreInit = TRUE)
  observeEvent(data.nuc2(), {
    x(sample(which(is.na(data.nuc2()$Class)),1))
  })
  observeEvent(input$previous, {
    x(x_list()[[length(x_list())]])
    x_list(x_list()[1:(length(x_list())-1)])
  })
  observeEvent(input$crop_n, {
    x(as.numeric(input$crop_n))
  })

  img.n = reactiveVal(NULL)
  obj.n = reactiveVal(NULL)
  class = reactiveVal(NULL)
  observeEvent(x(), {
    img.n(data.nuc2()$ImageNumber[x()])
    obj.n(data.nuc2()$ObjectNumber[x()])
  })
  
  output$x = renderText(x())
  output$img_n = renderText(img.n())
  output$obj_n = renderText(obj.n())
  output$img_name = renderText(unlist(data.nuc2()$Filename2[x()]))
  
  observeEvent(input$im2_dblclick, {
    brush(input$im2_brush)
  })
  output$im2 = renderImage({
    req(input$exp, x())
    img = image_read(paste0("MultiChannel_Images/Netosis_", input$exp, "_Merged_", data.nuc2()$Filename2[x()], ".jpg"))
    img2 = image_draw(img)
    coord = data.nuc2()[x(), c("X", "Y")]
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
    req(x(), obj.n(), img.n())
    im_list = map(list("Green" = "Green", "Phase" = "Phase", "Red" = "Red"), function(channel){
      filename = paste(data.nuc2()$Experiment[x()], channel, data.nuc2()$Filename2[x()], sep = "_")
      im = image_read(paste0("Raw_Images/", filename, ".png"))
      image_crop(im, geometry_area(50, 50, data.nuc2()$X[x()], data.nuc2()$Y[x()]))
    })
    im_Phase = image_level(im_list$Phase, black_point = input$Phase[1], white_point = input$Phase[2], mid_point = input$Phase2)
    im_Red = image_level(im_list$Red, black_point = input$Red[1], white_point = input$Red[2], mid_point = input$Red2)
    im_Green = image_level(im_list$Green, black_point = input$Green[1], white_point = input$Green[2], mid_point = input$Green2)
    (image_ggplot(im_Phase) + ggtitle("Phase")) / 
      (image_ggplot(im_Red) + ggtitle("Red")) / 
      (image_ggplot(im_Green) + ggtitle("Green")) & theme(plot.title = element_text(size = 20, hjust = 0.5))
  })
  output$n_class1 = renderTable({
    req(data.nuc2())
    table(data.nuc2()$Class) %>%
      as.data.frame %>% `colnames<-`(c("Class","Curated crops"))
  }, striped = TRUE)
  
  ## Tab Classified
  output$img_class = renderPlot({
    n.class = sum(data.nuc2()$Class == input$class, na.rm = TRUE)
    idx.class = sample(which(data.nuc2()$Class == input$class), min(60, n.class))
    crop.list = map(as.list(idx.class), function(z){
      img.z = image_read(paste0("MultiChannel_Images/Netosis_", input$exp, "_Merged_", data.nuc2()$Filename2[z], ".jpg"))
      crop.z = image_crop(img.z, geometry_area(width = 50, height = 50, x_off = data.nuc2()$X[z], y_off = data.nuc2()$Y[z]))
      image_ggplot(crop.z) +
        annotate("label", x = 43, y = 3, label = z, size = 3, hjust = 0.5)
    })
    wrap_plots(crop.list, ncol = 10, nrow = 6)
  })
  output$n_class2 = renderTable({
    req(data.nuc2())
    table(data.nuc2()$Class) %>%
      as.data.frame %>% `colnames<-`(c("Class","Curated crops"))
  }, striped = TRUE)
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))

