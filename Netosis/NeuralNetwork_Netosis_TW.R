## NN analyses cell crops
source("C:/Users/mwildschut/OneDrive - Vifor Pharma AG/Documents/R/Projects/R_CSL-Vifor_TW/SourceFile_TW.R")
setwd("C:/Incucyte/Netosis/Netosis_Exp10/Classification")
data = fread("Netosis_Exp10_Classification_2023-08-07_09-54-14.csv") %>%
  filter(!is.na(Class)) %>%
  group_by(Class) %>% add_tally(name = "Cells") %>% dplyr::slice_sample(n = min(.$Cells[.$Class != "NoCell"]))
map(as.list(unique(data$Class)), function(class){
  data.crops = data %>% filter(Class == class)
  for(i in 1:nrow(data.crops){
    for(channel in c("Phase", "Green", "Red"){
      
    })
    file.copy(paste0("Cell_Crops/Netosis-10_Phase_", data.nuc2()$Filename2[x()], "_", obj.n(), ".tif"))[1]
              im2 = image_read(paste0("Cell_Crops/Netosis-10_Red_", data.nuc2()$Filename2[x()], "_", obj.n(), ".tif"))[1]
              im3 = image_read(paste0("Cell_Crops/Netosis-10_Green_", data.nuc2()$Filename2[x()], "_", obj.n(), ".tif"))[1])
  }
})

# reticulate::install_miniconda()
# tensorflow::install_tensorflow()
keras::install_keras()
library(keras)
library(reticulate)
library(tensorflow)
tensorflow::
reticulate::use_miniconda(miniconda_path())
use_condaenv("r-tensorflow")
imdb <- dataset_imdb()

np <- import("numpy") #<- dataset_mnist()
mnist = np$load("mnist.npz")
mnist$files

a = py_to_r(mnist)
a = mnist$files
crop.size = 50

## Preparing the data 
c(c(x_train, y_train), c(x_test, y_test)) %<-% dataset_boston_housing()
c(c(train_images, train_labels), c(test_images, test_labels)) %<-% dataset_mnist()

dataset_mnist(path = "C:/Incucyte/Netosis/Netosis_Exp10/Classificationmnist.npz")
mnist$f[["data"]]


## Building the model --------------------------------------------------------------------------------------------------------------
model <- keras_model_sequential()%>%
  # Start with a hidden 2D convolutional layer
  layer_conv_2d(
    filters = 16, kernel_size = c(3,3), padding = "same",
    input_shape = c(crop.size, crop.size, 3), activation = 'leaky_relu'
  ) %>%
  
  # 2nd hidden layer
  layer_conv_2d(filter = crop.size, kernel_size = c(3,3), activation = 'leaky_relu') %>%
  
  
  # Use max pooling
  layer_max_pooling_2d(pool_size = c(2,2)) %>%
  layer_dropout(0.25) %>%
  
  # 3rd and 4th hidden 2D convolutional layers
  layer_conv_2d(filter = crop.size, kernel_size = c(3,3), padding = "same", activation = 'leaky_relu') %>%
  
  layer_conv_2d(filter = crop.size*2, kernel_size = c(3,3), activation = 'leaky_relu') %>%
  
  # Use max pooling
  layer_max_pooling_2d(pool_size = c(2,2)) %>%
  layer_dropout(0.25) %>%
  
  # Flatten max filtered output into feature vector
  # and feed into dense layer
  layer_flatten() %>%
  layer_dense(256, activation = 'leaky_relu') %>%
  layer_dropout(0.5) %>%
  
  # Outputs from dense layer
  layer_dense(6, activation = 'softmax')


## Compiling the model ----------------------------------------------------------------------------------------------------
learning_rate <- learning_rate_schedule_exponential_decay(
  initial_learning_rate = 5e-3,
  decay_rate = 0.96,
  decay_steps = 1500,
  staircase = TRUE
)
opt <- optimizer_adamax(learning_rate = learning_rate)

loss <- loss_sparse_categorical_crossentropy(from_logits = TRUE)

model %>% compile(
  loss = loss,
  optimizer = opt,
  metrics = "accuracy"
)


## Training the model ----------------------------------------------------------------------------------------------------
history <- model %>% fit(
  x_train, y_train,
  batch_size = 32,
  epochs = 10,
  validation_data = list(x_test, y_test),
  shuffle = TRUE
)


## Evaluating the model --------------------------------------------------------------------------------------------------
model %>% evaluate(x_test, y_test)
plot(history)
