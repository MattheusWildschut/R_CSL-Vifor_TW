## Read & pre-process images from Incucyte for CellSens/ScanR NNs

library(magick)

folder = "C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/Incucyte/Netosis/Netosis_Exp12/"
filename = "Netosis_Exp12_Green_A11_1_00d09h30m.tif"
files = list.files(paste0(folder, "Raw_Images"))
suppressWarnings(dir.create(paste0(folder, "Preprocessed_Images")))

walk(as.list(files), function(filename){
  exp = str_extract(filename, ".*(?=_Green|_Red|_Phase)")
  channel = str_extract(filename, "Green|Red|Phase")
  well = str_extract(filename, "[A-Z][:digit:]{2}")
  field = str_extract(filename, "(?<=[A-Z][:digit:]{2}_)[:digit:]*")
  mat96 = expand.grid(LETTERS[1:8], 1:12)
  mat96 = matrix(paste0(mat96$Var1, mat96$Var2), nrow = 8, ncol = 12)
  well2 = which(well == mat96)
  time = as.numeric(str_extract(filename, "(?<=d)[:digit:]{2}(?=h)"))*2 + ifelse(str_detect(filename, "30m"), 1, 0)
  fill.n = function(str, n = 5) {paste0(paste(rep(0, 5-nchar(str)), collapse = ""), str)}
  
  file.new = paste0(exp, "_", well, "--", "W", fill.n(well2), "_", "P", fill.n(field), "_", "T", fill.n(time), "--", channel, ".tif")
  
  im = image_read(paste0(folder, "Raw_Images/", filename))[1]
  im = image_scale(im, "x1000")
  image_write(im, path = paste0(folder, "Preprocessed_Images/", file.new), format = "tiff")
  print(paste0(file.new, " written to disk (", which(filename == files), "/", length(files), ")"))
})


filename = "Netosis_Exp12_Green_A11_1_00d09h30m.tif"
filename2 = "Netosis_Exp12_Red_A11_1_00d09h30m.tif"
filename3 = "Netosis_Exp12_Phase_A11_1_00d09h30m.tif"

im = image_read(paste0(folder, "Raw_Images/", filename))
im2 = image_read(paste0(folder, "Raw_Images/", filename2))
im3 = image_read(paste0(folder, "Raw_Images/", filename3))

im.c = image_join(im[1], im2[1], im3[1])
image_write(im.c, path = paste0(folder, "Preprocessed_Images/", "Merged_A11_8.tif"), 
            format = "tiff", depth = 16, channel = c("Red", "Green", "Grey"))

im = image_write(im, path = paste0(folder, "Preprocessed_Images/", file.new), format = "tiff")




