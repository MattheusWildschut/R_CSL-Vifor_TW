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
# in.folder = "E:/Netosis_Exp12/Renamed_Images"
in.folder = "C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/Incucyte/Netosis/Netosis_Exp10/Renamed_Images"
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
         Time = str_extract(Filename, "[:digit:]{2}d[:digit:]{2}h[:digit:]{2}m"),
         Time2 = as.numeric(factor(Time)))
files.new = paste0(f$Exp, "_", f$Well, "--", "W", fill.n(f$Well2), "--", "P", fill.n(f$Field), "--Z00000--", 
                   "T", fill.n(f$Time2), "--", f$Channel, ".tif")
# files.new = str_remove(files, "_[:digit:]{2}d[:digit:]{2}h[:digit:]{2}m")

## To rename -----------------------------------------------------------------------------------------------------------------------
rename = file.rename(paste0(in.folder, "/", files), paste0(in.folder, "/", files.new))
paste0(sum(rename), "/", length(files.new))

## To rewrite as 1024px ------------------------------------------------------------------------------------------------------
out.folder = "E:/Netosis_Exp12/Raw_Images_Renamed"
unlink(out.folder, recursive = TRUE)
dir.create(out.folder)

walk2(as.list(files), as.list(files.new), function(x,y){
  im = image_read(paste0(in.folder, "/", x), strip = TRUE)[1]
  im = image_scale(im, "x1024")
  image_write(im, path = paste0(out.folder, "/", y), format = "tiff", depth = 16)
  print(paste0(y, " written to disk (", which(x == files), "/", length(files), ")"))
})


## Trial multichannel TIFF -------------------------------------------------------------------------------------
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

## Analyze & plot Netosis data -----------------------------------------------------------------------------
folder = "C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/Incucyte/Netosis/Netosis_Exp12/Tables_ScanR"
data.class = fread(paste0(folder, "/All_NET-NN_Images-Converted4_NET_NN.txt")) %>%
  arrange(Area)
data.main = fread(paste0(folder, "/All_NET-NN_Images-Converted4_Main.txt"))
conditions = fread("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/Incucyte/Netosis/Netosis_Exp12/Conditions_Netosis-Exp12.csv") %>%
  mutate(Condition = str_remove(paste(Stimuli, Inhibitor, sep = "|"), "\\|$|^\\|")) %>%
  arrange(rank(as.numeric((str_extract(Inhibitor, "[:digit:]*"))), na.last = FALSE)) %>%
  arrange(rank(as.numeric((str_extract(Stimuli, "[:digit:]*"))), na.last = FALSE))
cond = unique(conditions$Condition)

mat96 = expand.grid(LETTERS[1:8], 1:12)
mat96 = matrix(paste0(mat96$Var1, mat96$Var2), nrow = 8, ncol = 12)
classes = c("NET", "SmallGreen", "Original", "FlatCell")
timepoints = paste0(rep(0:9, each = 2), c("h00m", "h30m"))

data.merged = data.main %>%
  left_join(data.class, by = join_by("Object ID" == "Parent Object ID (MO)"), multiple = "first") %>%
  # filter(Time>2) %>%
  mutate(Well = factor(t(mat96)[`Parent Object ID (Well)`+1], levels = as.vector(t(mat96))),
         Time2 = timepoints[Time]) %>%
  rowwise() %>% 
  dplyr::mutate(Max_prob = max(c(`Mean Intensity NET`, `Mean Intensity SmallGreen`, `Mean Intensity Original`, `Mean Intensity FlatCells`)),
                Class = ifelse(Max_prob == 0, 0, 
                               which(Max_prob == c(`Mean Intensity NET`, `Mean Intensity SmallGreen`, `Mean Intensity Original`, `Mean Intensity FlatCells`))),
                Class2 = factor(ifelse(Class == 0, "No class", classes[Class]), levels = c(classes, "No class"))) %>%
  filter(!is.na(Class2))
data.sum = data.merged %>%
  mutate(Condition = factor(conditions$Condition[match(Well, conditions$Well)], levels = cond),
         Class3 = factor(Class2, levels = c(classes, "No class"), labels = c("NET", "NET", "Original", "Original", "No class"))) %>%
  group_by(Condition, Well, Time2, Class2) %>% dplyr::summarise(N = n()) %>%
  group_by(Condition, Well, Time2) %>% dplyr::mutate(Percentage = round(N/sum(N)*100,1))# %>%
  # filter(!str_detect(Condition, "GSK"))

data.sum3 = data.sum %>%
  group_by(Condition, Time2, Class2) %>% dplyr::summarise(Percentage = round(mean(Percentage),1), N = sum(N))
# ggplotly(
a = as.character(unique(data.sum3$Condition))
             
  ggplot(data.sum3, aes(x = Condition, y = Percentage, fill = Class2)) +
    geom_col(position = position_stack(), col = "white", lwd = 0.1) +
    geom_vline(xintercept = c(7, 13), linetype = 2) +
    facet_wrap(.~Time2, dir = "h", ncol = 5) +
    scale_y_continuous(expand = expansion(mult = c(0,0)), labels = ~paste0(.x, "%")) +
    scale_fill_manual(values = c("NET" = "green3", "SmallGreen" = "green4", "Original" = "red4", "FlatCell" = "red1", "No class" = "grey")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), panel.grid = element_blank()) +
    labs(fill = "Class", x = NULL) +
    scale_x_discrete(limits = c(a[1:6], "", a[7:11], "", a[12:16]))
  # )

celltypes = c("NET", "SmallGreen")#c("FlatCell") #c("Original")#, "FlatCell") #c("NET", "SmallGreen")
data.sum2 = filter(data.sum, Class2 %in% celltypes)
# ggplotly(
plot.list = map(list("Medium", "PMA", "Ionomycin"), function(x){
  if(x == "Medium"){
    data.plot = data.sum %>% filter(!str_detect(Condition, "PMA|Ionomycin"))
  } else {
    data.plot = data.sum %>% filter(str_detect(Condition, x))
  }
  data.plot = data.plot %>% 
    filter(Class2 %in% celltypes & Condition != "20µM CMP") %>% 
    dplyr::mutate(Title = x, 
                  Condition = factor(Condition, levels = levels(Condition), 
                                     labels = str_replace(str_remove(levels(Condition), "5µM Ionomycin\\||100nM PMA\\|"),
                                                          "5µM Ionomycin|100nM PMA", "Medium")))
  ggplot(data.plot, aes(x = Time2, y = Percentage, col = Condition, group = Condition)) +
    # stat_summary(fun.data = mean_sdl, geom = "errorbar") + #fun.args = list(mult = 1),
    stat_summary(fun.y = ~sum(.x)/length(.x)*length(celltypes), geom = "line") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                       plot.title = element_text(hjust = 0.5, colour = "grey")) +
    facet_grid(.~Title) +
    scale_y_continuous(limits = c(0,85), expand = expansion(mult = c(0,0)), labels = ~paste0(.x, "%")) +
    labs(x = NULL, y = ifelse(x == "Medium", "Percentage NETs", ""))
})
wrap_plots(plot.list, nrow = 1, guides = "collect")  
# )

data.plot = data.sum %>%
  filter(Class2 %in% celltypes & Condition != "20µM CMP") %>%
  dplyr::mutate(Stimuli = factor(Condition, levels = levels(Condition),
                                   labels = replace_na(str_extract(levels(Condition), "PMA|Ionomycin"), "No stimuli")),
                Inhibitor = factor(Condition, levels = levels(Condition),
                                 labels = replace_na(str_extract(levels(Condition), "JBI_589|GSK484"), "GSK484")),
                Concentration = as.numeric(replace_na(str_extract(Condition, "62.5|125|250|500"), "0")))
ggplot(data.plot, aes(x = Time2, y = Percentage, col = Concentration, group = Concentration)) +
  # stat_summary(fun.data = mean_sdl, geom = "errorbar") + #fun.args = list(mult = 1),
  stat_summary(fun.y = ~sum(.x)/length(.x)*length(celltypes), geom = "line") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                     plot.title = element_text(hjust = 0.5, colour = "grey"), panel.grid = element_blank()) +
  facet_grid(Stimuli~Inhibitor) +
  scale_y_continuous(limits = c(0,82), expand = expansion(mult = c(0,0)), labels = ~paste0(.x, "%")) +
  scale_color_viridis(option = "H") +
  labs(x = NULL, "Percentage NETs")

## Analysis Exp10 -------------------------------------------------------------------------------------------------------------------
folder = "C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/Incucyte/Netosis/Netosis_Exp10/Tables_ScanR"
data.class = fread(paste0(folder, "/Netosis_Exp10_NET_NN.txt")) %>%
  arrange(Area)
data.main = fread(paste0(folder, "/Netosis_Exp10_Main.txt"))
# conditions = fread("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/Incucyte/Netosis/Netosis_Exp12/Conditions_Netosis-Exp12.csv") %>%
#   mutate(Condition = str_remove(paste(Stimuli, Inhibitor, sep = "|"), "\\|$|^\\|")) %>%
#   arrange(rank(as.numeric((str_extract(Inhibitor, "[:digit:]*"))), na.last = FALSE)) %>%
#   arrange(rank(as.numeric((str_extract(Stimuli, "[:digit:]*"))), na.last = FALSE))
conditions = fread(paste0(folder, "/PlateMap_Exp10.csv")) %>%
  mutate(Condition2 = paste(Condition, Concentration, sep = "_")) %>%
  arrange(Condition, Concentration)
cond = unique(conditions$Condition2)

mat96 = expand.grid(LETTERS[1:8], 1:12)
mat96 = matrix(paste0(mat96$Var1, mat96$Var2), nrow = 8, ncol = 12)
classes = c("NET", "SmallGreen", "Original", "FlatCell")
timepoints = paste0(rep(0:9, each = 2), c("h00m", "h30m"))

data.merged = data.main %>%
  left_join(data.class, by = join_by("Object ID" == "Parent Object ID (MO)"), multiple = "first") %>%
  # filter(Time>2) %>%
  mutate(Well = factor(t(mat96)[`Parent Object ID (Well)`+1], levels = as.vector(t(mat96))),
         Time2 = timepoints[Time]) %>%
  rowwise() %>% 
  dplyr::mutate(Max_prob = max(c(`Mean Intensity NET`, `Mean Intensity SmallGreen`, `Mean Intensity Original`, `Mean Intensity FlatCells`)),
                Class = ifelse(Max_prob == 0, 0, 
                               which(Max_prob == c(`Mean Intensity NET`, `Mean Intensity SmallGreen`, `Mean Intensity Original`, `Mean Intensity FlatCells`))),
                Class2 = factor(ifelse(Class == 0, "No class", classes[Class]), levels = c(classes, "No class"))) %>%
  filter(!is.na(Class2))
data.sum = data.merged %>%
  mutate(Condition = factor(conditions$Condition2[match(Well, conditions$Well)], levels = cond),
         Class3 = factor(Class2, levels = c(classes, "No class"), labels = c("NET", "NET", "Original", "Original", "No class"))) %>%
  group_by(Condition, Well, Time2, Class2) %>% dplyr::summarise(N = n()) %>%
  group_by(Condition, Well, Time2) %>% dplyr::mutate(Percentage = round(N/sum(N)*100,1))# %>%