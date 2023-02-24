## Data loading
setwd("C:\\Users\\mwildschut\\OneDrive - Vifor Pharma AG\\Documents\\R\\Projects\\R_CSL-Vifor_TW")
source("SourceFile_TW.R")

library(GEOquery)
library(DT)
dataset = "GSE7869"
gse = getGEO(dataset, GSEMatrix = TRUE)[[1]]

annot = gse@featureData@data
data = gse@assayData$exprs[annot$`Gene Symbol` != "",]
annot = annot[annot$`Gene Symbol` != "",]
annot.pat = gse@phenoData@data %>%
  mutate(Group = factor(case_when(str_detect(title, "normal|Normal") ~ "Normal",
                                  str_detect(title, "minimally") ~ "MCD",
                                  str_detect(title, "small") ~ "Cysts-Small",
                                  str_detect(title, "medium") ~ "Cysts-Medium",
                                  str_detect(title, "large") ~ "Cysts-Large"),
                        levels = c("Normal", "MCD", "Cysts-Small", "Cysts-Medium", "Cysts-Large")))
annot.genes = data.frame("Gene" = sort(unique(annot$`Gene Symbol`)), 
                         "Description" = annot$`Gene Title`[match(sort(unique(annot$`Gene Symbol`)), annot$`Gene Symbol`)])

## UI ----------------------------------------------------------------------------------------------------------------------------------
ui = {fluidPage(
  titlePanel("CSL-Vifor_Transcriptomics-browser_TW"),
  tabsetPanel(
    ## Data grouping UI --------
    tabPanel("Data grouping", fluid = TRUE, 
             sidebarLayout(
               sidebarPanel(
                 fluidRow(
                   radioGroupButtons(inputId = "grouping",
                                     label = "Select automatic or manual grouping",
                                     choices = c("Automatic", "Manual"),
                                     selected = "Automatic",
                                     justified = TRUE,
                                     checkIcon = list(yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                                                      no = tags$i(class = "fa fa-circle-o", style = "color: steelblue"))),
                   uiOutput(outputId = "ui.grouping"),
                   uiOutput(outputId = "group.name"),
                   column(6, uiOutput(outputId = "group.button")),
                   column(6, uiOutput(outputId = "reset.button")))
               ),
               mainPanel(
                 dataTableOutput("groups")
               )
             )
    ),
    ## Single gene analysis UI --------
    tabPanel("Single gene analysis", fluid = TRUE, 
             sidebarLayout(
               sidebarPanel(
                 selectizeInput(inputId = "gene",
                                label = "Select gene of interest",
                                choices = NULL),
                 materialSwitch(inputId = "probe.sum",
                                label = "Summarize probes per gene (max)",
                                value = TRUE,
                                status = "primary"),
                 uiOutput(outputId = "ui.probe"),
                 downloadBttn(outputId = "download.boxplot", 
                              label = list(icon("file-pdf"), "Download boxplot as PDF"),
                              color = "danger", size = "sm"),
                 hr(style = "border-color: black; color: black;"),
                 radioGroupButtons(inputId = "correlations",
                                   label = "Show correlations to selected gene",
                                   choices = c("Positive", "Negative", "Manual"),
                                   selected = "Positive",
                                   justified = TRUE,
                                   checkIcon = list(yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                                                    no = tags$i(class = "fa fa-circle-o", style = "color: steelblue"))),
                 uiOutput(outputId = "ui.manualcor"),
                 uiOutput(outputId = "ui.manualprobe"),
                 downloadBttn(outputId = "download.correlations", 
                              label = list(icon("file-pdf"), "Download correlations as PDF"),
                              color = "danger", size = "sm"),
                 hr(style = "border-color: black; color: black;"),
                 tags$b("Analyze correlations with Gene Set Enrichment Analysis (GSEA)"),
                 actionBttn(inputId = "GSEA",
                            label = "Run GSEA (takes a while)",
                            style = "unite", color = "success", size = "sm"),
                 downloadBttn(outputId = "download.corrGSEA", 
                              label = list(icon("file-pdf"), "Download GSEA as PDF"),
                              color = "danger", size = "sm")
               ),
               mainPanel(
                 fluidRow(
                   column(5, plotOutput("probe.boxplot")),
                   column(7, tags$b("Gene description:"), textOutput("gene.description"), HTML("<br/>"), 
                          DT::dataTableOutput("gene.ontologies")),
                   column(12, plotOutput("probe.correlations")),
                   column(12, plotOutput("probe.corrGSEA"))
                 )
               )
             )
    ), selected = "Single gene analysis" #"Data grouping" 
  )
)}

## Server ----------------------------------------------------------------------------------------------------------------------------------
server = function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })
  ## Data grouping Server ----------------------------------------------------
  output$ui.grouping = renderUI({
    if(input$grouping == "Manual"){
        pickerInput(
          inputId = "sample.grouping",
          label = "Select samples for group", 
          choices = annot.pat$geo_accession,
          choicesOpt = list(subtext = annot.pat$title),
          multiple = TRUE,
          options = list(
            `actions-box` = TRUE))
    }
  })
  output$group.name = renderUI({
    if(input$grouping == "Manual"){
      textInput(inputId = "group.name",
              label = "Write group name")
    }
  })
  output$group.button = renderUI({
    req(input$group.name)
    if(input$grouping == "Manual" & input$group.name != ""){
        actionBttn(inputId = "add.group",
               label = paste("Add group", input$group.name),
               style = "unite", 
               color = "primary")
    }
  })
  output$reset.button = renderUI({
    if(input$grouping == "Manual" & !is.null(unlist(reactiveValuesToList(group_list)))){
      actionBttn(inputId = "reset.group",
                 label = "Reset grouping",
                 style = "unite", 
                 color = "danger")
    }
  })
  group_list = reactiveValues()
  observeEvent(input$grouping, {
    map(names(group_list), function(i) group_list[[i]] = NULL)
    if(input$grouping == "Automatic"){
      map(as.list(levels(annot.pat$Group)), function(i) group_list[[i]] = annot.pat$geo_accession[annot.pat$Group == i])
    }
  })
  observeEvent(input$add.group, {
    group_list[[input$group.name]] <- c(group_list$dList, input$sample.grouping)
    updateTextInput(session, "group.name", value = "")
    updatePickerInput(session, "sample.grouping", selected = "")
  })
  observeEvent(input$reset.group, {
    map(names(group_list), function(i) group_list[[i]] = NULL)
    updateTextInput(session, "group.name", value = "")
    updatePickerInput(session, "sample.grouping", selected = "")
  })
  output$groups = DT::renderDataTable({
    req(input$grouping)
    if(!is.null(unlist(reactiveValuesToList(group_list)))){
      datatable({
        tab.groups = plyr::ldply(reactiveValuesToList(group_list), rbind) %>% `colnames<-`(c("Group name", paste("Sample", 1:(ncol(.)-1))))
        if(input$grouping == "Automatic") {
          tab.groups = tab.groups[match(levels(annot.pat$Group), tab.groups$`Group name`),]
        }
        tab.groups}, 
        rownames = FALSE, options = list(searching = FALSE, paging = FALSE, info = FALSE)) %>%
        formatStyle("Group name", fontWeight = "bold")
    }
  })
  
  ## Single gene analysis Server ----------------------------------------------------------------------------
  updateSelectizeInput(session, "gene", choices = annot.genes$Gene, selected = "SLC9B2", options = list(maxOptions = 100000), server = TRUE)
  output$gene.description = renderText({
    annot$`Gene Title`[annot$`Gene Symbol` == input$gene][1]
  })
  output$gene.ontologies = DT::renderDataTable({
    map(list("Biological Process", "Cellular Component", "Molecular Function"), function (x) {
      unlist(str_split(annot[,paste("Gene Ontology", x)][annot$`Gene Symbol` == input$gene][1], " /// ")) %>%
        str_split_fixed(" // ", 3) %>% as.data.frame %>% select(2) %>%
        t %>% as.data.frame
    }) %>% data.table::rbindlist(fill = TRUE) %>% t %>% as.data.frame %>% 
      `colnames<-`(c("Biological Process", "Cellular Component", "Molecular Function"))
  }, rownames = FALSE, options = list(pageLength = 6))
  output$ui.probe = renderUI({
    if (length(annot$ID[annot$`Gene Symbol` == input$gene])>1 & !input$probe.sum){
      selectInput(inputId = "select.probe",
                  label = paste0("Multiple probes found for gene ", input$gene, ": Select probe of interest"),
                  choices = annot$ID[annot$`Gene Symbol` == input$gene])}
  })
  select_probe = reactive({
    req(input$gene)
    if(!input$probe.sum) req(input$select.probe)
    if(input$probe.sum) {input$gene} else{
      ifelse(length(annot$ID[annot$`Gene Symbol` == input$gene])>1, input$select.probe, annot$ID[annot$`Gene Symbol` == input$gene])
    }
  })
  dataset = reactive({
    if(!input$probe.sum) return(data)
    as.data.frame(data) %>%
      mutate(Gene = annot$`Gene Symbol`) %>%
      group_by(Gene) %>%
      dplyr::summarise(across(everything(), max)) %>%
      column_to_rownames("Gene") %>% as.data.frame
  })
  annot.group = reactive({
    imap(reactiveValuesToList(group_list), function(x,y){
      temp = ifelse(match(colnames(dataset()), x)>0, y, "")
      temp[is.na(temp)] = ""
      temp
    }) %>% purrr::reduce(paste0)
  })
  probe_boxplot = reactive({
    req(input$gene)
    if (length(annot$ID[annot$`Gene Symbol` == input$gene])>1 & !input$probe.sum){
      req(input$select.probe)
    }
    data.probe = as.data.frame(t(dataset()[select_probe(),,drop = FALSE])) %>%
      mutate(Group = if(input$grouping == "Automatic") annot.pat$Group else annot.group()) %>%
      filter(Group != "") %>%
      `colnames<-`(c("Probe", "Group"))
    ggplot(data.probe, aes(x = Group, y = Probe, fill = Group)) +
      geom_boxplot() +
      labs(x = NULL, y = paste0(input$gene, " (", select_probe(), ")\nlog2 mRNA levels")) +
      theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  })
  output$ui.manualcor = renderUI({
    req(input$gene)
    if(input$correlations == "Manual"){
      selectizeInput(inputId = "corrgene",
                     label = paste0("Select gene of interest to correlate against ", input$gene),#, " (", select_probe(), ")"),
                     choices = NULL)
      }
    })
  output$ui.manualprobe = renderUI({
    if (input$correlations == "Manual" & length(annot$ID[annot$`Gene Symbol` == input$corrgene])>1 & !input$probe.sum){
      selectInput(inputId = "select.corrprobe",
                  label = paste0("Multiple probes found for gene ", input$corrgene, ": Select probe of interest"),
                  choices = annot$ID[annot$`Gene Symbol` == input$corrgene])
    }
  })
  observe({
    if(input$correlations == "Manual"){
      updateSelectizeInput(session = session, inputId = "corrgene", choices = annot.genes$Gene, selected = "PDGFRB",
                           options = list(maxOptions = 100000), server = TRUE)
    }
  })
  select_corrprobe = reactive({
    if(input$correlations != "Manual") return(NULL)
    if(input$probe.sum) return(input$corrgene)
    req(input$corrgene, input$select.corrprobe)
    ifelse(length(annot$ID[annot$`Gene Symbol` == input$corrgene])>1,
           input$select.corrprobe, annot$ID[annot$`Gene Symbol` == input$corrgene])
  })
  cors = reactive({
    req(input$gene, select_probe())
    if (length(annot$ID[annot$`Gene Symbol` == input$gene])>1 & !input$probe.sum){
      req(input$select.probe)
    }
    group = if(input$grouping == "Automatic") annot.pat$Group else annot.group()
    if(input$probe.sum){
      annot = annot %>% group_by(`Gene Symbol`) %>% dplyr::summarise(ID = unique(`Gene Symbol`)) %>%
        mutate(`Gene Title` = annot$`Gene Title`[match(ID, annot$`Gene Symbol`)])
    } 
    as.data.frame(cor(t(dataset()[,group != ""]), unlist(dataset()[select_probe(), group != ""]))) %>%
      `colnames<-`("Correlation") %>%
      mutate(ID = rownames(.),
             Gene = annot$`Gene Symbol`[match(annot$ID, ID)],
             Description = annot$`Gene Title`[match(annot$ID, ID)]) %>%
      filter(rownames(.) != select_probe())
  })
  probe_correlations = reactive({
    if(input$correlations == "Manual" & !input$probe.sum) req(input$select.corrprobe, input$corrgene)
    if(input$correlations == "Manual") req(input$corrgene)
      cor_function = function(direction, i){
      cors.sel = if(direction == "Negative") cors() %>% arrange(Correlation) else cors() %>% arrange(-Correlation)
      idx = if(direction == "Manual") filter(cors.sel, ID == select_corrprobe()) else cors.sel[i,]
      data.cor = cbind.data.frame(unlist(dataset()[select_probe(),]), unlist(dataset()[idx$ID,])) %>%
        `colnames<-`(c("Input", "Gene")) %>%
        mutate(Group = if(input$grouping == "Automatic") annot.pat$Group else annot.group()) %>%
        filter(Group != "")
      scatter.plot = ggplot(data.cor, aes(x = Input, y = Gene)) +
        geom_point(aes(col = Group), show.legend = FALSE) +
        geom_smooth(method = "lm", formula = y ~ x) + stat_cor() +
        labs(x = ifelse(input$probe.sum, paste(input$gene, "(Max probe)"), paste0(input$gene, " (", select_probe(), ")")), 
             y = ifelse(input$probe.sum, paste(idx$Gene, "(Max probe)"), paste0(idx$Gene, " (", idx$ID, ")"))) +
        theme_bw()
      box.plot = ggplot(data.cor, aes(x = Group, y = Gene)) +
        geom_boxplot(aes(fill = Group)) +
        labs(x = NULL, y = ifelse(input$probe.sum, paste(idx$Gene, "(Max probe)"), paste0(idx$Gene, " (", idx$ID, ")"))) +
        theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
      return(list("scatter.plot" = scatter.plot, "box.plot" = box.plot))
    }
    n.plots = if(input$correlations == "Manual") 1 else 1:5
    scatter.list = map(n.plots, function(x) cor_function(input$correlations, x)$scatter.plot)
    boxplot.list = map(n.plots, function(x) cor_function(input$correlations, x)$box.plot)
    wrap_plots(scatter.list, nrow = 1) / wrap_plots(boxplot.list, nrow = 1) + plot_layout(guides = "collect")
  })
  probe_corrGSEA = eventReactive(input$GSEA, {
    if(input$probe.sum){
      annot = annot %>% group_by(`Gene Symbol`) %>% dplyr::summarise(ID = unique(`Gene Symbol`)) %>%
        mutate(`Gene Title` = annot$`Gene Title`[match(ID, annot$`Gene Symbol`)],
               ENTREZ_GENE_ID = annot$ENTREZ_GENE_ID[match(ID, annot$`Gene Symbol`)])
    } 
    cors.sum = cors() %>%
      mutate(ENTREZID = str_remove(annot$ENTREZ_GENE_ID[match(ID, annot$ID)], " ///.*")) %>%
      filter(ENTREZID != "") %>%
      group_by(ENTREZID) %>%
      dplyr::summarise(Correlation = mean(Correlation)) %>%
      arrange(-Correlation)
    
    cors.list = setNames(cors.sum$Correlation, cors.sum$ENTREZID)
    set.seed(0)
    enr.KEGG = gseKEGG(cors.list, keyType = "ncbi-geneid", organism = "hsa", eps = 0)
    data.plot = rbind.data.frame(slice_max(enr.KEGG@result, NES, n = 10), slice_min(enr.KEGG@result, NES, n = 10)) %>%
      mutate(Description = factor(Description, levels = Description),
             Sign = factor(sign(NES), levels = c(-1,1), labels = c("Negative", "Positive")))

    ggplot(data.plot, aes(x = NES, y = Description, size = -log10(`p.adjust`), col = setSize)) +
      geom_point() + scale_color_viridis() +
      facet_wrap(vars(Sign), scales = "free") +
      labs(x = "NES (Normalized enrichment score)", y = NULL,
           title = paste0("Top 10 enriched KEGG terms\nfrom gene set enrichment analysis (GSEA)\non ", 
                          input$gene, " (", ifelse(input$probe.sum, "Max probe", select_probe()), ") correlations")) +
      scale_y_discrete(limits = rev) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  })
  output$probe.boxplot = renderPlot({
    probe_boxplot()
  })
  output$probe.correlations = renderPlot({
    probe_correlations()
  }, width = reactive(ifelse(input$correlations == "Manual", 325, "auto")))
  output$probe.corrGSEA = renderPlot({
    probe_corrGSEA()
  })
  output$download.boxplot = downloadHandler(
    filename = function() { paste0('Boxplot_', input$gene, "-", select_probe(), "_", Sys.Date(), '.pdf') },
    content = function(file) {
      ggsave(file, probe_boxplot(), height = 6, width = 5)
  })
  output$download.correlations = downloadHandler(
    filename = function() { paste0('Correlations_', input$gene, "-", select_probe(), "_", Sys.Date(), '.pdf') },
    content = function(file) {
      ggsave(file, probe_correlations(), height = 6, width = ifelse(input$correlations == "Manual", 4.75, 16))
  })
  output$download.corrGSEA = downloadHandler(
    filename = function() { paste0('Correlation-GSEA_', input$gene, "-", select_probe(), "_", Sys.Date(), '.pdf') },
    content = function(file) {
      ggsave(file, probe_corrGSEA(), height = 6, width = 12)
  })
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))
               