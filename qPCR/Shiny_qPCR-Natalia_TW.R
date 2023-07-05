setwd("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/qPCR")
source("../SourceFile_TW.R")
library(shinyFiles)

## UI ----------------------------------------------------------------------------------------------------------------------------------
ui = fluidPage(
  # useShinyjs(),
  titlePanel("Natalia-qPCR_TW"),
  sidebarLayout(
    sidebarPanel(
      shinyFilesButton("data_files", label = "Select data file(s)", title = "Please select data file(s)", multiple = TRUE,
                       buttonType = "success"),
      verbatimTextOutput("data_filepaths"),
      
      shinyFilesButton("layout_file", label = "Select layout file", title = "Please select layout file", multiple = FALSE,
                       buttonType = "primary"),
      # verbatimTextOutput('rawInputValue'),
      verbatimTextOutput("layout_filepath"),
      uiOutput("ui.layout")
      #,
      #verbatimTextOutput('rawInputValue')
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data table", dataTableOutput("data"),
                 uiOutput("ui.download.xlsx")),
        tabPanel("Plot | Influence medium", uiOutput("medium.ui"),
                 # put sliders & download button as UI that only shows up when plots is ready?
                 sliderInput("medium.columns", label = "Select amount of columns", min = 1, max = 10, value = 4),
                 sliderInput("medium.width", label = "Select width plot", min = 100, max = 900, value = 625, step = 25), 
                 sliderInput("medium.height", label = "Select height plot", min = 100, max = 900, value = 425, step = 25)),
        tabPanel("Plot | qPCR data", uiOutput("plot.ui"),
                 sliderInput("plot.columns", label = "Select amount of columns", min = 1, max = 10, value = 4),
                 sliderInput("plot.width", label = "Select width plot", min = 100, max = 1500, value = 1250, step = 25), 
                 sliderInput("plot.height", label = "Select height plot", min = 100, max = 900, value = 500, step = 25))
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
  
  # output$rawInputDataValue <- renderPrint({str(input$data_files)})
  output$data_filepaths = renderText({
    if(is.integer(input$data_files)){
      "Please select data file(s)"
    } else {
      paste0("Selected data files:\n",
             paste(paste0(1:length(data_filepaths()), ". ", data_filepaths()), collapse = "\n"))
    }
  })
  # output$rawInputValue <- renderPrint({str(input$layout_file)})
  output$layout_filepath = renderText({
    if(is.integer(input$layout_file)){
      "Please select layout file"
    } else {
      paste0("Selected layout file:\n", layout_filepath())
    }
  })
  output$data = renderDataTable({
    req(!is.integer(input$data_files))
    map_dfr(as.list(data_filepaths()), fread)
  })
  output$layout = renderDataTable({
    req(!is.integer(input$layout_file))
    data.frame(readxl::excel_sheets(path = layout_filepath()))
  })
  output$ui.layout = renderUI({
    req(!is.integer(input$layout_file))
    fluidRow(column(radioGroupButtons(inputId = "layout_sheet",
                                      label = "Select layout sheet",
                                      choices = readxl::excel_sheets(path = layout_filepath()),
                                      size = "sm", individual = TRUE, selected = character(0),
                                      checkIcon = list(yes = icon("circle", class = "fa-solid fa-circle", style = "color: steelblue"),
                                                       no = icon("circle", class = "fa-regular fa-circle", style = "color: steelblue"))), 
                    width = 12),
             column(radioGroupButtons(inputId = "sample_sheet",
                                      label = "Select sample list sheet",
                                      choices = readxl::excel_sheets(path = layout_filepath()),
                                      size = "sm", individual = TRUE, selected = character(0),
                                      checkIcon = list(yes = icon("circle", class = "fa-solid fa-circle", style = "color: steelblue"),
                                                       no = icon("circle", class = "fa-regular fa-circle", style = "color: steelblue"))), 
                    width = 12)
    )
  })
  # output$rawInputValue <- renderPrint({str(input$sample_sheet)})
  
  data = reactive({
    req(!is.integer(input$data_files) & !is.integer(input$layout_file) & !is.null(input$layout_sheet) & !is.null(input$sample_sheet))
    
    raw.data1 = fread(data_filepaths()[1])
    raw.data2 = fread(data_filepaths()[2])

    sample.list = readxl::read_excel(layout_filepath(), sheet = input$sample_sheet)
    raw.layout = xlsx::read.xlsx(layout_filepath(), sheetName = input$layout_sheet, header = FALSE)

    layout = raw.layout %>%
      `colnames<-`(LETTERS[1:ncol(.)]) %>%
      dplyr::mutate(across(everything(), as.character)) %>%
      rownames_to_column("Row") %>% pivot_longer(-Row, names_to = "Column", values_to = "Condition") %>%
      dplyr::mutate(Cell = paste0(Column, Row))
    layout.genes = raw.layout %>% dplyr::select(1) %>% unlist %>% .[!is.na(.) & . != "Gene"]

    { ## Function to get cell colors out of Excel file
      wb     <- loadWorkbook(layout_filepath())
      sheet1 <- getSheets(wb)[[input$layout_sheet]]
      styles <- sapply(getCells(getRows(sheet1)), getCellStyle)

      cellColor <- function(style)
      {
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

      colors.genes = colors.cells %>% filter(Column.num == 1 & Row != 1) %>%
        mutate(Gene = layout.genes) %>%
        dplyr::select(Cell, Gene, Color)
    }

    layout.function = function(rows, columns){
      colors.cells %>% filter(Row %in% rows & str_detect(Cell, columns)) %>%
        mutate(Gene = colors.genes$Gene[match(Color, colors.genes$Color)],
               Condition.num = layout$Condition[match(Cell, layout$Cell)],
               Condition = factor(sample.list$Samples[match(Condition.num, sample.list$No)], levels = sample.list$Samples),
               Column.num = Column.num-min(Column.num)+1,
               Row = Row-min(Row)+1,
               Well = paste0(LETTERS[Row], Column.num))
    }

    data1 = layout.function(rows = 2:17, columns = "[C-Z]") %>%
      left_join(raw.data1 %>% select(Pos, Cp), by = join_by("Well" == "Pos"))
    data2 = layout.function(rows = 21:36, columns = "[C-Z]") %>%
      left_join(raw.data2 %>% select(Pos, Cp), by = join_by("Well" == "Pos"))

    data = rbind.data.frame(data1, data2) %>%
      filter(!is.na(Condition) & !Condition %in% c("NTC", "Vehicle")) %>%
      group_by(Condition) %>% dplyr::mutate(Cp_norm = round(Cp-mean(Cp[Gene == "HPRT"], na.rm = TRUE),2)) %>%
      filter(Gene != "HPRT") %>% mutate(Gene = factor(Gene, levels = layout.genes)) %>% arrange(Gene) %>%
      group_by(Gene) %>% dplyr::mutate(Gene2 = factor(Gene, labels = unique(paste0(Gene, " (Ct: ", round(mean(Cp[Condition == "Medium"]),1), ")")))) %>%
      group_by(Gene) %>% dplyr::mutate(mRNA_norm = 2^-(Cp_norm-mean(Cp_norm[Condition == "Medium"], na.rm = TRUE))) %>%
      mutate(Medium = factor(replace_na(str_extract(Condition, "\\+ Cytotox Green"), "Medium"),
                             levels = c("Medium", "+ Cytotox Green")), #labels = c("Medium", "Medium + Cytotox Green")),
             Condition2 = factor(str_remove(Condition, " \\(\\+ Cytotox Green\\)"),
                                 levels = unique(str_remove(sample.list$Samples, " \\(\\+ Cytotox Green\\)")))) %>%
      group_by(Gene, Medium) %>% dplyr::mutate(`Normalized mRNA expression` = round(2^-(Cp_norm-mean(Cp_norm[str_detect(Condition, "Medium")], na.rm = TRUE)),2)) %>%
      arrange(Gene, Condition, `Normalized mRNA expression`)
  })
  output$data = renderDataTable({
    data() %>% dplyr::select(Cell, Well, Gene, Condition, Medium, Cp, Cp_norm, `Normalized mRNA expression`)
  })
  output$medium = renderPlot({
    data.medium = data() %>% filter(Condition2 == "Medium")
    ggplot(data.medium, aes(x = Condition, y = mRNA_norm, fill = Condition, shape = Medium)) +
      geom_hline(yintercept = 1, linetype = 2) +
      geom_point(size = 2.5, col = "grey30", stroke = 0.01) +
      facet_wrap(.~Gene, ncol = input$medium.columns) +
      scale_shape_manual(values = c("Medium" = 21, "+ Cytotox Green" = 23)) +
      theme_bw() + theme(text = element_text(size = 18), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                         panel.grid.major.x = element_blank()) +
      labs(x = NULL, y = "mRNA expression\n(normalized to medium)") + #scale_y_continuous(trans = "log2") + #limits = c(0, NA), expand = expansion(mult = c(0.03,0.05)), 
      guides(fill = guide_legend(override.aes = list(shape = c(21,23))),
             shape = "none")
  })
  output$medium.ui = renderUI({
    plotOutput("medium", width = input$medium.width, height = input$medium.height)
  })
  output$plot = renderPlot({
    req(input$sample_sheet)
    sample.list = readxl::read_excel(layout_filepath(), sheet = input$sample_sheet)
    ggplot(data(), aes(x = Condition, y = `Normalized mRNA expression`, fill = Condition2, shape = Medium)) +
      geom_hline(yintercept = 1, linetype = 2) +
      geom_point(size = 2.5, col = "grey30", stroke = 0.01) +
      geom_vline(xintercept = nrow(sample.list)/2+0.5, linetype = 2) + 
      facet_wrap(.~Gene2, scales = "free", ncol = input$plot.columns) +
      scale_shape_manual(values = c("Medium" = 21, "+ Cytotox Green" = 23)) +
      scale_fill_discrete(type = c("black", brewer.pal("Paired", n = length(unique(data()$Condition2))-1))) +
      theme_bw() + theme(text = element_text(size = 18), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                         panel.grid.major.x = element_blank()) +
      labs(x = NULL, fill = "Condition") + scale_y_continuous(trans = "log2") + #limits = c(0, NA), expand = expansion(mult = c(0.03,0.05)), 
      guides(fill = guide_legend(override.aes = list(shape = 21)),
             shape = guide_legend(override.aes = list(fill = "black")))
  })
  output$plot.ui = renderUI({
    plotOutput("plot", width = input$plot.width, height = input$plot.height)
  })
  output$ui.download.xlsx = renderUI({
    req(!is.integer(input$data_files) & !is.integer(input$layout_file), input$layout_sheet, input$sample_sheet)
    downloadBttn(outputId = "download.xlsx",
                 label = list(icon("file-excel"), "Download data as .xlsx"),
                 color = "success", size = "sm")
  })
  output$download.xlsx = downloadHandler(
    filename = function() { paste0('Excel-data_qPCR_', Sys.Date(), '.xlsx') },
    content = function(file) {
      wb = createWorkbook()
      addDataFrame(rbind(c("Folder", getwd()),
                         c("Data file(s)", paste(data_filepaths(), collapse = " | ")),
                         c("Layout file", layout_filepath()),
                         c("Layout sheet", input$layout_sheet),
                         c("Sample sheet", input$sample_sheet)),
                   sheet = createSheet(wb, "Metadata"), row.names = FALSE, col.names = FALSE)
      addDataFrame(as.data.frame(data() %>% dplyr::select(Cell, Well, Gene, Condition, Medium, Cp, Cp_norm, `Normalized mRNA expression`)), 
                   sheet = createSheet(wb, "qPCR_Data"), row.names = FALSE)
      saveWorkbook(wb, file)
  })
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))