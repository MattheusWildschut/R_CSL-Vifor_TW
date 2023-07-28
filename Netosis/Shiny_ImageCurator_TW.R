## Shiny Image Curator
setwd("C:/Incucyte/Netosis/")
source("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/SourceFile_TW.R")
library(shinyFiles)
library(magick)

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
         Y = ifelse(Location_Center_Y-(crop.size/2) > 0, Location_Center_Y-(crop.size/2), 0),
         Border = X < 0 | X > 1264-crop.size | Y < 0 | Y > 936-crop.size,
         Class = as.character(NA)) %>%
  filter(!Border)

## UI ----------------------------------------------------------------------------------------------------------------------------------
ui = fluidPage(
  titlePanel("Shiny Image Curator_TW"),
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("data_folder", label = "Select folder", title = "Please select folder", buttonType = "success"),
      verbatimTextOutput("folder_path")
    ),
    mainPanel(
      column(width = 6, plotOutput("im")),
      column(width = 6, imageOutput("im2")),
      dataTableOutput("tab"),
      textOutput("value"),
      textOutput("img.n"),
      textOutput("obj.n"),
      column(width = 4,
             sliderInput(inputId = "Phase", label = "Phase min/max threshold", min = 0, max = 100, value = c(0,100), step = 1),
             sliderInput(inputId = "Phase2", label = "Phase midpoint", min = 0, max = 2, value = 1, step = 0.01)),
      column(width = 4,
             sliderInput(inputId = "Red", label = "Red min/max threshold", min = 0, max = 2, value = c(0,0.75), step = 0.01),
             sliderInput(inputId = "Red2", label = "Red midpoint", min = 0, max = 2, value = 0.5, step = 0.01)),
      column(width = 4,
             sliderInput(inputId = "Green", label = "Green min/max threshold", min = 0, max = 25, value = c(0,2), step = 0.01),
             sliderInput(inputId = "Green2", label = "Green midpoint", min = 0, max = 2, value = 0.5, step = 0.01)),
      # actionBttn(inputId = "Next", label = "Next image", style = "unite", color = "succes"),
      # actionBttn(inputId = "Previous", label = "Previous image", style = "unite", color = "danger"),
      actionBttn(inputId = "NET", label = "NET", style = "unite", color = "succes")
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
  # data.nuc = reactiveVal({
  #   data.nuc
  # })
  data.nuc2 = reactiveVal(data.nuc)
  x = reactiveVal(1)
  img.n = reactiveVal(1)
  obj.n = reactiveVal(1)
  observeEvent(input$NET, {
    data.nuc = data.nuc2()
    data.nuc2(data.nuc %>%
                mutate(Class = ifelse(ObjectNumber == obj.n() & ImageNumber == img.n(), "NET", Class)))
    x(sample(which(is.na(data.nuc2()$Class)),1))
    img.n(data.nuc$ImageNumber[x()])
    obj.n(data.nuc$ObjectNumber[x()])
  })
  output$img.n = renderText(img.n())
  output$obj.n = renderText(obj.n())
  output$tab = renderDataTable({
    data.nuc2()
  })
  # x = reactiveVal(1)
  # observeEvent(input$Next, {
  #   x(x()+1)
  # })
  # observeEvent(input$Previous, {
  #   x(x()-1)
  # })
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
    # img = image_read(paste0("C:/Incucyte/Netosis/Netosis_Exp10/MultiChannel_Images/Netosis_Exp10_Merged_", data.nuc$Filename2[x()], ".jpg"))
    img = image_read(paste0("C:/Incucyte/Netosis/Netosis_Exp10/MultiChannel_Images/Netosis_Exp10_Merged_", "C10_1_00d08h00m", ".jpg"))
    img2 = image_draw(img)
    coord = data.nuc[x(), c("X", "Y")]
    rect(coord$X, coord$Y, coord$X+50, coord$Y+50, border = "white", lty = "dashed", lwd = 3)
    outfile <- tempfile(fileext='.png')
    image_write(img2, outfile)
    list(src = outfile, width = 540,  height = 400)
  }, deleteFile = TRUE)
}

## App ----------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server, options = list("launch.browser" = TRUE))

img = image_read("C:/Incucyte/Netosis/Netosis_Exp10/Incucyte_Images/Netosis_Exp10_Merged_C10_1_00d08h00m.jpg")
img2 = image_draw(img)
coord = data.nuc %>% filter(ObjectNumber == 1 & ImageNumber == 1) %>% select(X, Y)
rect(coord$X, coord$Y, coord$X+50, coord$Y+50, border = "white", lty = "dashed", lwd = 3)
dev.off()



# im2 = image_read(paste0(folder, "/Netosis-10_Red_C10_1_00d08h00m_", 1, ".tif"))[1]
# image_normalize(im)
# im2 = image_level(im2, black_point = 0, white_point = .75, mid_point = .5, channel = "Red")
# image_ggplot(im2, interpolate = TRUE)
# image_level(im2, black_point = 0, white_point = 2, mid_point = .5, channel = "Red")
# 
# image_colorize(im2, 20, "red")
# image_fill(im2, "red")
# 
# a = image_ggplot(im2)
# image_ggplot(im2) + scale_fill_viridis()
# a$scales$scales
# 
# im2j = im2 %>% image_write(tempfile(fileext='.jpg'), format = 'jpg')
# image_read(im2j, channel = "red")
# 
# image_draw(im2, pointsize = 2, res = .5)
# a = as.data.frame(as.integer(im2[[1]]))
# col_fun = colorRamp2(c(min(a), max(a)), c("black", "red"))
# red = as_ggplot(grid.grabExpr(draw(Heatmap(as.matrix(as.data.frame(as.integer(im2[[1]]))), col = col_fun, 
#                                            cluster_rows = FALSE, cluster_columns = FALSE,
#                                            show_column_names = FALSE, show_heatmap_legend = FALSE))))
# 
# as_ggplot(grid.grabExpr(draw(Heatmap(matrix(1:10, ncol =2)))))
# dev.off()
# dev.new()
# 
# library(EBImage)
# img = readImage(paste0(folder, "/Netosis-10_Red_C10_1_00d08h00m_", 1, ".tif"))
# plot_ly(x = 1:50, y = 1:50, z = img@.Data, type="heatmap", colors = colorRamp(c("black", "red")))
# img = readImage('https://upload.wikimedia.org/wikipedia/commons/thumb/0/00/Crab_Nebula.jpg/240px-Crab_Nebula.jpg')
# volcano
# 
# library(ijtiff)
# img <- array(seq_len(2^3), dim = c(2, 2, 2, 1)) #  2 channel, 1 frame
# ijtiff_img(img, description = "blah blah")
# 
# library(raster)
# folder = "C:/Incucyte/Netosis/Netosis_Exp10/Cell_Crops"
# red = brick(raster(paste0(folder, "/Netosis-10_Red_C10_1_00d08h00m_", 1, ".tif")))
# grayscale_colors <- gray.colors(100,            # number of different color levels 
#                                 start = 0.0,    # how black (0) to go
#                                 end = 1.0,      # how white (1) to go
#                                 gamma = 2.2,    # correction between how a digital 
#                                 # camera sees the world and how human eyes see it
#                                 alpha = NULL)   #Null=colors are not transparent
# plot(red, col = grayscale_colors, axes = FALSE)
# green = brick(raster(paste0(folder, "/Netosis-10_Green_C10_1_00d08h00m_", 1, ".tif")))
# phase = brick(raster(paste0(folder, "/Netosis-10_Phase_C10_1_00d08h00m_", 1, ".tif")))
# plot(red)
# all = stack(red, green, phase)
# plotRGB(all, scale = 50000)#, stretch = "lin")
# 
# draster = raster(as(cells,"SpatialPixelsDataFrame"))
# b = brick(draster, 1-draster, 1-draster)
# plotRGB(b, scale=1)
# plot(trig, col=NA, border='white', lwd=5, add=T)
# 
# 
# im <- abs(matrix(rnorm(100), 10, 10))
# mI <- max(im, na.rm=TRUE)
# rM <- im 
# rM[is.na(rM)] <- 0
# colRamp <- colorRamp(c("black", "red"))
# col <- rgb(colRamp(rM/mI), maxColorValue=255)
# 
# p <- matrix(NA, nrow=nrow(im), ncol=ncol(im), byrow=TRUE)
# 
# p[,] <- col
# # 
# # layout(matrix(c(1,2), nrow=2, ncol=1), heights=c(4,1))
# # layout.show(2)
# # par(mar=c(1,1,1,1))
# # breaks <- seq(min(rM/mI), max(rM/mI),length.out=100)
# # 
# plot(NA, type="n", xlim=c(1, ncol(p)), ylim=c(1, nrow(p)),
#      xlab="", ylab="", xaxs="i", yaxs="i", axes=FALSE, asp=1)
# rasterImage(as.raster(p), xleft=1, xright=ncol(p), ybottom=1, ytop=nrow(p),
#             interpolate=TRUE)
# 
# par(mar=c(3,1,1,1)) 
# pal.1=colorRampPalette(c("black", "blue", "green", "yellow", "red", "purple"), space="rgb")
# image.scale(m3, col=pal.1(length(breaks)-1), breaks=breaks, horiz=TRUE)
# 
# img = image_read("C:/Incucyte/Netosis/Netosis_Exp10/Incucyte_Images/Netosis_Exp10_Merged_C10_1_00d08h00m.jpg")
# img2 = image_draw(img)
# rect(20, 50, 70, 100, border = "white", lty = "dashed", lwd = 3)
# dev.off()
# print(img2)
