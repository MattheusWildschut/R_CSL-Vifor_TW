folder = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(folder)
Sys.setenv(LANG = "EN")
suppressPackageStartupMessages(suppressWarnings({
  library(stringr)
  library(viridis)
  library(ggplot2)
  library(data.table)
  library(ggpubr)
  library(shiny)
  library(shinyWidgets)
  library(shinyFiles)
  library(knitr)
  library(patchwork)
  library(ggfortify)
  library(xlsx)
  library(tidyverse)
  library(shinyjs)
  library(sortable)
}))
options(dplyr.summarise.inform = FALSE)

## UI ----------------------------------------------------------------------------------------------------------------------------------
ui = fluidPage(
  titlePanel("Natalia-qPCR_TW"),
  sidebarLayout(
    sidebarPanel(
      shinyFilesButton("data_files", label = "Select data file(s)", title = "Please select data file(s)", multiple = TRUE,
                       buttonType = "success"),
      verbatimTextOutput("data_filepaths"),
      shinyFilesButton("layout_file", label = "Select layout file", title = "Please select layout file", multiple = FALSE,
                       buttonType = "primary"),
      verbatimTextOutput("layout_filepath"),
      uiOutput("ui.layout"),
      uiOutput("ui.layout2")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data table", dataTableOutput("data.table"),
                 uiOutput("ui.download.xlsx")),
        tabPanel("Plot | Influence medium", uiOutput("medium.ui"), uiOutput("medium.ui2")),
        tabPanel("Plot | qPCR data", uiOutput("plot.ui"), uiOutput("plot.ui2")),
        tabPanel("Plot | Delta CT values", uiOutput("dCTplot.ui"), uiOutput("dCTplot.ui2"))
      )
    )
  )
)

## Server ----------------------------------------------------------------------------------------------------------------------------------
server = function(input, output, session) {
  session$onSessionEnded(function() {stopApp()})
  roots = c("Input" = "Input/", Home = "../")
  shinyFileChoose(input, "data_files", roots=roots, filetypes=c("", "txt", "csv", "xls", "xlsx"))
  data_filepaths <- reactive(unlist(parseFilePaths(roots, input$data_files)[,"datapath"]))
  shinyFileChoose(input, "layout_file", roots=roots, filetypes=c("", "txt", "csv", "xls", "xlsx"))
  layout_filepath <- reactive(unlist(parseFilePaths(roots, input$layout_file)[,"datapath"]))
  
  output$data_filepaths = renderText({
    if(is.integer(input$data_files)){
      "Please select data file(s)"
    } else {
      paste0("Selected data files:\n",
             paste(paste0(1:length(data_filepaths()), ". ", data_filepaths()), collapse = "\n"))
    }
  })
  output$layout_filepath = renderText({
    if(is.integer(input$layout_file)){
      "Please select layout file"
    } else {
      paste0("Selected layout file:\n", layout_filepath())
    }
  })
  layout_sheets = reactive(as.list(1:length(data_filepaths())))
  output$ui.layout = renderUI({
    req(!is.integer(input$layout_file))
    fluidRow(column(radioGroupButtons(inputId = "sample_sheet",
                                      label = "Select sample list sheet",
                                      choices = readxl::excel_sheets(path = layout_filepath())[str_detect(readxl::excel_sheets(path = layout_filepath()),"(?i)sample")],
                                      size = "sm", individual = TRUE, selected = character(0),
                                      checkIcon = list(yes = icon("circle", class = "fa-solid fa-circle", style = "color: steelblue"),
                                                       no = icon("circle", class = "fa-regular fa-circle", style = "color: steelblue"))), 
                    width = 12),
             map(layout_sheets(), function(i){
               column(radioGroupButtons(inputId = paste0("layout_sheet", i),
                                      label = paste0("Select layout sheet for ", "'", data_filepaths()[i], "'"),
                                      choices = readxl::excel_sheets(path = layout_filepath())[str_detect(readxl::excel_sheets(path = layout_filepath()),"(?i)layout")],
                                      size = "sm", individual = TRUE, selected = character(0),
                                      checkIcon = list(yes = icon("circle", class = "fa-solid fa-circle", style = "color: steelblue"),
                                                       no = icon("circle", class = "fa-regular fa-circle", style = "color: steelblue"))), 
                    width = 12)})
             )
  })
  layout_sheets2 = reactive(map(as.list(1:length(data_filepaths())), ~input[[paste0("layout_sheet", .x)]]))
  data = reactive({
    req(!is.integer(input$data_files) & !is.integer(input$layout_file) & !is.null(input$sample_sheet) &
          all(unlist(map(layout_sheets2(), ~!is.null(.x)))))

    sample.list = readxl::read_excel(layout_filepath(), sheet = input$sample_sheet)
    
    data.comb = map2_dfr(as.list(data_filepaths()), as.list(layout_sheets2()), function(path, sheet){
      raw.layout = xlsx::read.xlsx(layout_filepath(), sheetName = sheet, header = FALSE) %>%
        `colnames<-`(paste0(rep(c("", LETTERS), each = length(LETTERS)), LETTERS)[1:ncol(.)])

      layout = raw.layout %>%
        dplyr::mutate(across(everything(), as.character)) %>%
        rownames_to_column("Row") %>% pivot_longer(-Row, names_to = "Column", values_to = "Condition") %>%
        dplyr::mutate(Cell = paste0(Column, Row))
      layout.genes = raw.layout %>% dplyr::select(1) %>% unlist %>% .[!is.na(.) & . != "Gene"]
  
      wb     <- loadWorkbook(layout_filepath())
      sheet1 <- getSheets(wb)[[sheet]]
      styles <- sapply(getCells(getRows(sheet1)), getCellStyle)

      cellColor <- function(style){
        fg  <- style$getFillForegroundXSSFColor()
        rgb <- tryCatch(fg$getRgb(), error = function(e) NULL)
        rgb <- paste(rgb, collapse = "")
        return(rgb)
      }
      colors.cells = sapply(styles, cellColor) %>% as.data.frame %>% `colnames<-`("Color") %>%
        rownames_to_column("Index") %>%
        dplyr::mutate(Column.num = as.numeric(str_extract(Index, "(?<=\\.)[:digit:]{1,2}")),
                      Column = LETTERS[Column.num],
                      Row = as.numeric(str_extract(Index, "[:digit:]{1,2}(?=\\.)")),
                      Cell = paste0(LETTERS[Column.num], Row)) %>%
        arrange(Column.num, Row) 

      colors.genes = colors.cells %>% filter(Column.num == 1 & Row != 1 & Color != "") %>%
        mutate(Gene = layout.genes) %>%
        dplyr::select(Cell, Gene, Color)
      
      colors.cells = colors.cells %>% 
        filter(Row %in% 2:17 & str_detect(Column, "[C-Z]")) %>%
        mutate(Gene = factor(colors.genes$Gene[match(Color, colors.genes$Color)], levels = layout.genes),
               Condition.num = layout$Condition[match(Cell, layout$Cell)],
               Condition = factor(sample.list$Samples[match(Condition.num, sample.list$No)],
                                  levels = c("Medium", unique(sample.list$Samples)[unique(sample.list$Samples) != "Medium"])),
               Medium = factor(sample.list$Medium[match(Condition.num, sample.list$No)],
                               levels = c("Medium", unique(sample.list$Medium[sample.list$Medium != "Medium"]))),
               Column.num = Column.num-min(Column.num)+1,
               Row = Row-min(Row)+1,
               Well = paste0(LETTERS[Row], Column.num))
      left_join(colors.cells, fread(path), by = join_by("Well" == "Pos")) %>%
        mutate(Well = paste0(str_extract(sheet, "Plate [:digit:]|Plate[:digit:]"), "_", Well))
    })
    
    data = data.comb %>%
      filter(!is.na(Condition) & !Condition %in% c("NTC", "Vehicle", "Incucyte - Vehicle")) %>%
      dplyr::mutate(Cp = ifelse(Cp == 40, NA, Cp)) %>%
      group_by(Condition, Medium) %>% dplyr::mutate(Cp_norm = round(Cp-mean(Cp[Gene == "HPRT"], na.rm = TRUE),2)) %>%
      arrange(Gene) %>%
      group_by(Gene) %>% dplyr::mutate(Gene2 = factor(Gene, labels = unique(paste0(Gene, " (Ct: ", round(mean(Cp[Condition == "Medium"], na.rm = TRUE),1), ")")))) %>%
      group_by(Gene) %>% dplyr::mutate(mRNA_norm = 2^-(Cp_norm-mean(Cp_norm[Condition == "Medium"], na.rm = TRUE))) %>%
      group_by(Gene, Medium) %>% dplyr::mutate(`Normalized mRNA expression` = round(2^-(Cp_norm-mean(Cp_norm[str_detect(Condition, "Medium")], na.rm = TRUE)),2)) %>%
      arrange(str_extract(Well, "[:alpha:]"), as.numeric(str_extract(Well, "[:digit:]+")))
  })
  output$data.table = renderDataTable({
    data() %>% dplyr::select(Well, Gene, Condition, Medium, Cp, Cp_norm, `Normalized mRNA expression`)
  })
  output$ui.layout2 = renderUI({
    req(data())
    colors = c("black", "steelblue1", "steelblue4", "firebrick1", "firebrick4", "springgreen1", "springgreen4", "darkorange1", "darkorange4", "pink1", "pink4", 
                      "steelblue2", "steelblue3", "firebrick2", "firebrick3", "springgreen2", "springgreen3", "darkorange2", "darkorange3", "pink2", "pink3",
                      "goldenrod", "violet", "magenta")
    fluidRow(tags$style(HTML(paste0("#", colors, " {background-color: ", col2hcl(colors), c("; color: white", rep("", length(colors)-1)), ";}"))),
             column(rank_list(input_id = "list_mediums",
                              text = tags$b("Sort medium conditions (first will be used as 'base')"),
                              labels = unique(data()$Medium)),
                    width = 12),
             column(rank_list(input_id = "list_conditions",
                              text = tags$b("Sort stimuli conditions"),
                              labels = unique(data()$Condition)),
                    width = 6),
             column(rank_list(input_id = "list_colors",
                              text = tags$b("Sort colors"),
                              labels = map(as.list(colors), function(x) div(x, id = x))),
                    width = 6))
  })
  medium_react = reactive({
    data.medium = data() %>% filter(Condition == "Medium") %>% 
      dplyr::mutate(Medium = factor(Medium, levels = input$list_mediums))
    ggplot(data.medium, aes(x = Medium, y = mRNA_norm, fill = Medium, shape = Medium)) +
      geom_hline(yintercept = 1, linetype = 2) +
      geom_point(size = 2.5, col = "grey30", stroke = 0.01) +
      facet_wrap(.~Gene, ncol = input$medium.columns) +
      scale_shape_manual(values = c(21, 23, 22, 24, 25, 26, 27)[1:length(levels(data.medium$Medium))]) +
      theme_bw() + theme(text = element_text(size = 18), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                         panel.grid.major.x = element_blank()) +
      labs(x = NULL, y = paste0("mRNA expression Medium condition\n(normalized to ", input$list_mediums[1], ")"))
  })
  output$medium = renderPlot({
    medium_react()
  })
  output$medium.ui2 = renderUI({
    req(!is.integer(input$data_files) & !is.integer(input$layout_file) &
          all(unlist(map(layout_sheets2(), ~!is.null(.x)))), input$sample_sheet)
    fluidRow(hr(), 
      column(sliderInput("medium.columns", label = "Select amount of columns", min = 1, max = 10, value = 4), hr(), 
             downloadBttn("medium.download", label = list(icon("file-pdf"), "Download plot as PDF"), color = "danger", size = "sm"), width = 3),
      column(sliderInput("medium.width", label = "Select width plot", min = 100, max = 900, value = 625, step = 25), hr(), width = 3),
      column(sliderInput("medium.height", label = "Select height plot", min = 100, max = 900, value = 425, step = 25), hr(), width = 3)
    )
  })
  output$medium.ui = renderUI({
    req(input$medium.width)
    plotOutput("medium", width = input$medium.width, height = input$medium.height)
  })
  plot_react = reactive({
    req(input$sample_sheet)
    sample.list = readxl::read_excel(layout_filepath(), sheet = input$sample_sheet)
    plot.data = data() %>% filter(Gene != "HPRT") %>%
      dplyr::mutate(Medium = factor(Medium, levels = input$list_mediums), Condition = factor(Condition, levels = input$list_conditions)) %>%
      arrange(Medium, Condition) %>% mutate(Xaxis = factor(paste(Medium, Condition), levels = unique(paste(Medium, Condition))))
    x_intercepts = floor(quantile(1:length(unique(plot.data$Xaxis)), seq(0,1,1/length(levels(plot.data$Medium)))[!seq(0,1,1/length(levels(plot.data$Medium))) %in% c(0,1)])) + 0.5
    ggplot(plot.data, aes(x = Xaxis, y = `Normalized mRNA expression`, fill = Condition, shape = Medium)) +
      geom_hline(yintercept = 1, linetype = 2) +
      geom_point(size = 2.5, col = "grey30", stroke = 0.01) +
      geom_vline(xintercept = x_intercepts, linetype = 2) + 
      facet_wrap(.~Gene2, scales = "free", ncol = input$plot.columns) +
      scale_shape_manual(values = c(21, 23, 22, 24, 25, 26, 27)[1:length(levels(plot.data$Medium))]) +
      scale_fill_discrete(type = input$list_colors[1:length(levels(plot.data$Condition))]) +
      theme_bw() + theme(text = element_text(size = 18), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                         panel.grid.major.x = element_blank()) +
      labs(x = NULL, fill = "Condition") + scale_y_continuous(trans = "log2") + 
      guides(fill = guide_legend(override.aes = list(shape = 21)),
             shape = guide_legend(override.aes = list(fill = "black")))
  })
  output$plot = renderPlot({
    plot_react()
  })
  output$plot.ui2 = renderUI({
    req(!is.integer(input$data_files) & !is.integer(input$layout_file) &
          all(unlist(map(layout_sheets2(), ~!is.null(.x)))), input$sample_sheet)
    fluidRow(hr(), 
      column(sliderInput("plot.columns", label = "Select amount of columns", min = 1, max = 10, value = 4), hr(), 
             downloadBttn("plot.download", label = list(icon("file-pdf"), "Download plot as PDF"), color = "danger", size = "sm"), width = 3),
      column(sliderInput("plot.width", label = "Select width plot", min = 150, max = 1500, value = 1250, step = 25), hr(), width = 3),
      column(sliderInput("plot.height", label = "Select height plot", min = 100, max = 900, value = 500, step = 25), hr(), width = 3)
    )
  })
  output$plot.ui = renderUI({
    req(input$plot.width)
    plotOutput("plot", width = input$plot.width, height = input$plot.height)
  })
  dCTplot_react = reactive({
    req(input$sample_sheet)
    sample.list = readxl::read_excel(layout_filepath(), sheet = input$sample_sheet)
    plot.data = data() %>% filter(Gene != "HPRT") %>%
      dplyr::mutate(Cp_norm = 2^-Cp_norm) %>%
      dplyr::mutate(Medium = factor(Medium, levels = input$list_mediums), Condition = factor(Condition, levels = input$list_conditions)) %>%
      arrange(Medium, Condition) %>% mutate(Xaxis = factor(paste(Medium, Condition), levels = unique(paste(Medium, Condition))))
    x_intercepts = floor(quantile(1:length(unique(plot.data$Xaxis)), seq(0,1,1/length(levels(plot.data$Medium)))[!seq(0,1,1/length(levels(plot.data$Medium))) %in% c(0,1)])) + 0.5
    # avg_medium = mean(plot.data$Cp_norm[plot.data$Medium == input$list_mediums[1] & plot.data$Condition == "Medium"], na.rm = TRUE)
    ggplot(plot.data, aes(x = Xaxis, y = Cp_norm, fill = Condition, shape = Medium)) +
      # geom_hline(yintercept = avg_medium, linetype = 2) +
      geom_point(size = 2.5, col = "grey30", stroke = 0.01) +
      geom_vline(xintercept = x_intercepts, linetype = 2) + 
      facet_wrap(.~Gene2, scales = "free", ncol = input$dCTplot.columns) +
      scale_shape_manual(values = c(21, 23, 22, 24, 25, 26, 27)[1:length(levels(plot.data$Medium))]) +
      scale_fill_discrete(type = input$list_colors[1:length(levels(plot.data$Condition))]) +
      theme_bw() + theme(text = element_text(size = 18), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                         panel.grid.major.x = element_blank()) +
      labs(x = NULL, fill = "Condition") + scale_y_continuous(trans = "log2") + 
      guides(fill = guide_legend(override.aes = list(shape = 21)),
             shape = guide_legend(override.aes = list(fill = "black")))
  })
  output$dCTplot = renderPlot({
    dCTplot_react()
  })
  output$dCTplot.ui2 = renderUI({
    req(!is.integer(input$data_files) & !is.integer(input$layout_file) &
          all(unlist(map(layout_sheets2(), ~!is.null(.x)))), input$sample_sheet)
    fluidRow(hr(), 
             column(sliderInput("dCTplot.columns", label = "Select amount of columns", min = 1, max = 10, value = 4), hr(), 
                    downloadBttn("dCTplot.download", label = list(icon("file-pdf"), "Download plot as PDF"), color = "danger", size = "sm"), width = 3),
             column(sliderInput("dCTplot.width", label = "Select width plot", min = 150, max = 1500, value = 1250, step = 25), hr(), width = 3),
             column(sliderInput("dCTplot.height", label = "Select height plot", min = 100, max = 900, value = 500, step = 25), hr(), width = 3)
    )
  })
  output$dCTplot.ui = renderUI({
    req(input$dCTplot.width)
    plotOutput("dCTplot", width = input$dCTplot.width, height = input$dCTplot.height)
  })
  output$ui.download.xlsx = renderUI({
    req(!is.integer(input$data_files) & !is.integer(input$layout_file) &
          all(unlist(map(layout_sheets2(), ~!is.null(.x)))), input$sample_sheet)
    downloadBttn(outputId = "download.xlsx",
                 label = list(icon("file-excel"), "Download data as .xlsx"),
                 color = "success", size = "sm")
  })
  output$download.xlsx = downloadHandler(
    filename = function() { paste0('Excel-data_qPCR_', Sys.Date(), '.xlsx') },
    content = function(file) {
      wb = createWorkbook()
      addDataFrame(rbind(c("Folder", getwd()),
                         c("Data file(s)", paste(data_filepaths(), collapse = " | ")), c("Layout file", layout_filepath()),
                         c("Layout sheet", paste(layout_sheets2(), collapse = " | ")), c("Sample sheet", input$sample_sheet)),
                   sheet = createSheet(wb, "Metadata"), row.names = FALSE, col.names = FALSE)
      addDataFrame(as.data.frame(data() %>% dplyr::select(Cell, Well, Gene, Condition, Medium, Cp, Cp_norm, `Normalized mRNA expression`)), 
                   sheet = createSheet(wb, "qPCR_Data"), row.names = FALSE)
      saveWorkbook(wb, file)
  })
  output$medium.download = downloadHandler(
    filename = function() { paste0("qPCR-Medium-plot_", Sys.Date(), '.pdf') },
    content = function(file) {
      ggsave(file, medium_react(), width = input$medium.width/72, height = input$medium.height/72)
  })
  output$plot.download = downloadHandler(
    filename = function() { paste0("qPCR-plot_", Sys.Date(), '.pdf') },
    content = function(file) {
      ggsave(file, plot_react(), width = input$plot.width/72, height = input$plot.height/72)
  })
  output$dCTplot.download = downloadHandler(
    filename = function() { paste0("qPCR-dCTplot_", Sys.Date(), '.pdf') },
    content = function(file) {
      ggsave(file, dCTplot_react(), width = input$dCTplot.width/72, height = input$dCTplot.height/72)
    })
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))
