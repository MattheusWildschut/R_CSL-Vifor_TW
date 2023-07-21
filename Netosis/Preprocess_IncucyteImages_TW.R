## Read & pre-process images from Incucyte for CellSens/ScanR NNs

source("SourceFile_TW.R")
library(magick)

# in.folder = "C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/Incucyte/Netosis/Netosis_Exp12/Raw_Images/"
in.folder = "E:/RPTEC3_Injury_Cytotox/Raw_Images"
files = list.files(in.folder)
# filename = "Netosis_Exp12_Green_A11_1_00d09h30m.tif"
# filename = files[175]

out.folder = "E:/RPTEC3_Injury_Cytotox/Raw_Images_Renamed1024px"
if(exists(out.folder)){
  unlink(out.folder, recursive = TRUE)
  dir.create(out.folder)
} else {dir.create(out.folder)}

fill.n = function(txt, n = 5) {paste0(paste(rep(0, 5-nchar(txt)), collapse = ""), txt)}

walk(as.list(files), function(filename){
  exp = str_extract(filename, ".*(?=_Green|_Red|_Phase)")
  channel = str_extract(filename, "Green|Red|Phase")
  well = str_extract(filename, "(?<=_)[A-Z][:digit:]{1,2}")
  field = str_extract(filename, "(?<=_[A-Z][:digit:]{1,2}_)[:digit:]*")
  mat96 = expand.grid(LETTERS[1:8], 1:12)
  mat96 = matrix(paste0(mat96$Var1, mat96$Var2), nrow = 8, ncol = 12)
  well2 = which(well == t(mat96))
  time = as.numeric(str_extract(filename, "(?<=d)[:digit:]{2}(?=h)"))*2 + ifelse(str_detect(filename, "30m"), 1, 0)
  
  file.new = paste0(exp, "_", well, "--", "W", fill.n(well2), "--", "P", fill.n(field), "--Z00000--", "T", fill.n(time), "--", channel, ".tif")
  
  im = image_read(paste0(in.folder, "/", filename), strip = TRUE)[1]
  im = image_scale(im, "x1024")
  image_write(im, path = paste0(out.folder, "/", file.new), format = "tiff", depth = 16)
  print(paste0(file.new, " written to disk (", which(filename == files), "/", length(files), ")"))
})

## New version ------------------------------------------------------------------------------------------
# in.folder = "C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/Incucyte/Netosis/Netosis_Exp12/Raw_Images_Renamed/"
in.folder = "F:/RPTEC3_Injury_Cytotox_Renamed/"
files = list.files(in.folder)
mat96 = expand.grid(LETTERS[1:8], 1:12)
mat96 = matrix(paste0(mat96$Var1, mat96$Var2), nrow = 8, ncol = 12)
fill.n = function(txt) unlist(map(as.list(txt), ~paste0(paste(rep(0, 5-nchar(.x)), collapse = ""), .x)))

f = files %>%
  as.data.frame %>% `colnames<-`("Filename") %>%
  mutate(Exp = str_extract(Filename, ".*(?=_Green|_Red|_Phase)"),
         Channel = str_extract(Filename, "Green|Red|Phase"),
         Well = str_extract(Filename, "(?<=_)[A-Z][:digit:]{1,2}"),
         Well2 = match(Well, t(mat96)),
         Field = str_extract(Filename, "(?<=_[A-Z][:digit:]{1,2}_)[:digit:]*"),
         Time = str_extract(Filename, "[:digit:]{2}d[:digit:]{2}h[:digit:]{2}m")) %>%
  group_by(Exp, Channel, Well, Field) %>% dplyr::mutate(Time2 = order(Time))
files.new = paste0(f$Exp, "_", f$Well, "--", "W", fill.n(f$Well2), "--", "P", fill.n(f$Field), "--Z00000--", 
                   "T", fill.n(f$Time2), "--", f$Channel, ".tif")

## To rename
rename = file.rename(paste0(in.folder, files), paste0(in.folder, files.new))
paste0(sum(rename), "/", length(files.new))

## To rewrite as 1024px
map2(as.list(files), as.list(files.new), function(x,y){
  im = image_read(paste0(in.folder, "/", x), strip = TRUE)[1]
  im = image_scale(im, "x1024")
  image_write(im, path = paste0(out.folder, "/", y), format = "tiff", depth = 16)
  print(paste0(y, " written to disk (", which(filename == files), "/", length(files), ")"))
})



TIFFWriteDirectoryTagData
magick:::magick_image_write()

files[c(175, 1020,1196,1535)]

filename = "Netosis_Exp12_Green_A11_1_00d09h30m.tif"
filename2 = "Netosis_Exp12_Red_A11_1_00d09h30m.tif"
filename3 = "Netosis_Exp12_Phase_A11_1_00d09h30m.tif"

im = image_read(paste0(folder, "Raw_Images/", filename))
im2 = image_read(paste0(folder, "Raw_Images/", filename2))
im3 = image_read(paste0(folder, "Raw_Images/", filename3))

im.c = image_join(im, im2, im3)
image_write(im.c, path = paste0(folder, "Preprocessed_Images/", "Merged_A11_2.tif"), 
            format = "tiff", depth = 16)

im = image_write(im, path = paste0(folder, "Preprocessed_Images/", file.new), format = "tiff")

image_data(im.c, channels = c("Red", "Green", "Grey"))

length(im.c)

image = im.c
channels = c("Red", "Green", "Grey")

img.function = function (image, channels = NULL, frame = 1) 
{
  if (length(image) > 1 || frame > 1) 
    image <- image[frame]
  if (!length(channels)) {
    info <- image_info(image)
    channels <- if (tolower(info$colorspace) == "gray") {
      "gray"
    }
    else if (isTRUE(info$matte)) {
      "rgba"
    }
    else {
      "rgb"
    }
  }
  image <- image_flatten(image)
  image_write_frame(image, format = channels)
}
img.function(im.c, c("Red", "Green", "Grey"))

im.f = image_flatten(im.c)

## 
data.class = fread(paste0(folder, "ParameterData_Classification.txt")) %>%
  arrange(Area)
data.main = fread(paste0(folder, "ParameterData_Main.txt"))

mat96 = expand.grid(LETTERS[1:8], 1:12)
mat96 = matrix(paste0(mat96$Var1, mat96$Var2), nrow = 8, ncol = 12)

data.merged = data.main %>%
  left_join(data.class, by = join_by("Object ID" == "Parent Object ID (MO)"), multiple = "first") %>%
  mutate(Well = t(mat96)[`Parent Object ID (Well)`+1],
         Class = c("NET", "SmallGreen", "Original", "FlatCell")[`Most frequent class Classification`]) %>%
  filter(!is.na(Class)) %>%
  group_by(Well, Class) %>% dplyr::summarise(N = n()) %>%
  group_by(Well) %>% dplyr::mutate(Frequency = N/sum(N))

ggplot(data.merged, aes(x = Class, y = Frequency, fill = Class)) +
  geom_col() +
  facet_wrap(.~Well, ncol = 3, dir = "v") +
  theme_bw() + scale_y_continuous(expand = expansion(mult = c(0,0.05)), labels = scales::percent) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + labs(x = NULL)


ggplotly(ggplot(data.img, aes(x = Timepoint, y = NET_percentage, col = Concentration, group = Concentration)) +
           stat_summary(fun.data = mean_sdl, geom = "errorbar") + #fun.args = list(mult = 1),
           stat_summary(fun.y = mean, geom = "line") +
           # geom_line() +
           facet_grid(.~Condition) +
           scale_color_viridis_d() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)))
