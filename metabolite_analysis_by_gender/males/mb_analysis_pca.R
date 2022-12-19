## IMPORTANT NOTE #######################################################################################
## Execute mb_data_filtering.R - From line 1 to 51 ######################################################
## To see PCA before batch effect correction, execute from line 1 to line 44 ############################
#########################################################################################################

# Load necessary libraries
pacman::p_load(caret, ggplot2)

attach(ds)

## Scale data - Use caret library to MinMax scale all columns
process <- preProcess(as.data.frame(ds), method = "range")
norm_scale <- predict(process, as.data.frame(ds))

## Batch-Effect PCA #####################################################################################
## Insert Batch column in ds dataframe and concatenate the scaled dataset
ds <- data.frame(dataset['Batch'], norm_scale)
## Change Batch to categorical variable
ds$Batch <- factor(ds$Batch, levels=c(0,1), labels=c("Batch1", "Batch2"))

## PCA experiment - Get dataset to work on, i.e. remove Batch which will be our variable
pca_dataset <- ds[ , 2:163]

## PCA ##################################################################################################
pca <- prcomp(pca_dataset)
pca

## returns the scree plot showing the variances explained by the principal ## components versus the principal component number
plot(pca, type = "l")

# Get summary of the pca variable
summary(pca)  # Cutoff after 3

## Plot PCA #############################################################################################
## 2D PCA ###############################################################################################
library(ggfortify)
pca_dataset$Batch <- ds$Batch
autoplot(pca, data = pca_dataset, colour = 'Batch', label = FALSE)

library(ggbiplot)
g1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = ds$Batch, ellipse = TRUE, circle = TRUE)
g1 <- g1 + scale_color_discrete(name = '')
g1 <- g1 + theme(legend.direction = 'horizontal', legend.position = 'top')
g1
#########################################################################################################

## Fractures PCA ########################################################################################
## Insert Status Fractures column in ds dataframe
ds <- data.frame(dataset['Status_Fractures'], norm_scale)
## Change Batch to categorical variable
ds$Status_Fractures <- factor(ds$Status_Fractures, levels=c(0,1), labels=c("No", "Yes"))

## PCA experiment - Get dataset to work on, i.e. remove Status_Fractures which will be our variable
pca_dataset <- ds[ , 2:163]

## PCA ##################################################################################################
pca <- prcomp(pca_dataset)
pca

## returns the scree plot showing the variances explained by the principal ## components versus the principal component number
plot(pca, type = "l")

# Get summary of the pca variable
summary(pca)  # Cutoff after 3

## Plot PCA #############################################################################################
## 2D PCA ###############################################################################################
library(ggfortify)
pca_dataset$Status_Fractures <- ds$Status_Fractures
autoplot(pca, data = pca_dataset, colour = 'Status_Fractures', label = FALSE)

library(ggbiplot)
g1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = ds$Status_Fractures, ellipse = TRUE, circle = TRUE)
g1 <- g1 + scale_color_discrete(name = '')
g1 <- g1 + theme(legend.direction = 'horizontal', legend.position = 'top')
g1
#########################################################################################################

## BMD PCA ##############################################################################################
## Insert Status BMD column in ds dataframe
ds <- data.frame(dataset['Status_BMD'], norm_scale)
## Change Status_BMD to categorical variable
ds$Status_BMD <- factor(ds$Status_BMD, levels=c(0,1,2), labels=c("Normal", "Osteopenic", "Osteoporotic"))

## PCA experiment - Get dataset to work on, i.e. remove Status_BMD which will be our variable
pca_dataset <- ds[ , 2:163]

## PCA ##################################################################################################
pca <- prcomp(pca_dataset)
pca

## returns the scree plot showing the variances explained by the principal ## components versus the principal component number
plot(pca, type = "l")

# Get summary of the pca variable
summary(pca)  # Cutoff after 3

## Plot PCA #############################################################################################
## 2D PCA ###############################################################################################
library(ggfortify)
pca_dataset$Status_BMD <- ds$Status_BMD
autoplot(pca, data = pca_dataset, colour = 'Status_BMD', label = FALSE)

library(ggbiplot)
g1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = ds$Status_BMD, ellipse = TRUE, circle = TRUE)
g1 <- g1 + scale_color_discrete(name = '')
g1 <- g1 + theme(legend.direction = 'horizontal', legend.position = 'top')
g1
#########################################################################################################

## Lumbar Spine PCA #####################################################################################
## Insert LS Status column in ds dataframe
ds <- data.frame(dataset['LS_status'], norm_scale)
## Change LS_Status to categorical variable
ds$LS_status <- factor(ds$LS_status, levels=c(0,1,2), labels=c("Normal", "Osteopenic", "Osteoporotic"))

## PCA experiment - Get dataset to work on, i.e. remove LS_Status which will be our variable
pca_dataset <- ds[ , 2:163]

## PCA ##################################################################################################
pca <- prcomp(pca_dataset)
pca

## returns the scree plot showing the variances explained by the principal ## components versus the principal component number
plot(pca, type = "l")

# Get summary of the pca variable
summary(pca)  # Cutoff after 3

## Plot PCA #############################################################################################
## 2D PCA ###############################################################################################
library(ggfortify)
pca_dataset$LS_status <- ds$LS_status
autoplot(pca, data = pca_dataset, colour = 'LS_status', label = FALSE)

library(ggbiplot)
g1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = ds$LS_status, ellipse = TRUE, circle = TRUE)
g1 <- g1 + scale_color_discrete(name = '')
g1 <- g1 + theme(legend.direction = 'horizontal', legend.position = 'top')
g1
#########################################################################################################

## Femoral Neck PCA #####################################################################################
## Insert FN Status column in ds dataframe
ds <- data.frame(dataset['FN_status'], norm_scale)
## Change FN_Status to categorical variable
ds$FN_status <- factor(ds$FN_status, levels=c(0,1,2), labels=c("Normal", "Osteopenic", "Osteoporotic"))

## PCA experiment - Get dataset to work on, i.e. remove FN_Status which will be our variable
pca_dataset <- ds[ , 2:163]

## PCA ##################################################################################################
pca <- prcomp(pca_dataset)
pca

## returns the scree plot showing the variances explained by the principal ## components versus the principal component number
plot(pca, type = "l")

# Get summary of the pca variable
summary(pca)  # Cutoff after 3

## Plot PCA #############################################################################################
## 2D PCA ###############################################################################################
library(ggfortify)
pca_dataset$FN_status <- ds$FN_status
autoplot(pca, data = pca_dataset, colour = 'FN_status', label = FALSE)

library(ggbiplot)
g1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = ds$FN_status, ellipse = TRUE, circle = TRUE)
g1 <- g1 + scale_color_discrete(name = '')
g1 <- g1 + theme(legend.direction = 'horizontal', legend.position = 'top')
g1
#########################################################################################################

## Total Hip PCA ########################################################################################
## Insert FN Status column in ds dataframe
ds <- data.frame(dataset['TH_status'], norm_scale)
## Change TH_Status to categorical variable
ds$TH_status <- factor(ds$TH_status, levels=c(0,1,2), labels=c("Normal", "Osteopenic", "Osteoporotic"))

## PCA experiment - Get dataset to work on, i.e. remove TH_Status which will be our variable
pca_dataset <- ds[ , 2:163]

## PCA ##################################################################################################
pca <- prcomp(pca_dataset)
pca

## returns the scree plot showing the variances explained by the principal ## components versus the principal component number
plot(pca, type = "l")

# Get summary of the pca variable
summary(pca)  # Cutoff after 3

## Plot PCA #############################################################################################
## 2D PCA ###############################################################################################
library(ggfortify)
pca_dataset$TH_status <- ds$TH_status
autoplot(pca, data = pca_dataset, colour = 'TH_status', label = FALSE)

library(ggbiplot)
g1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = ds$TH_status, ellipse = TRUE, circle = TRUE)
g1 <- g1 + scale_color_discrete(name = '')
g1 <- g1 + theme(legend.direction = 'horizontal', legend.position = 'top')
g1
#########################################################################################################


## Execute PCA on Fractures and BMD according to p-values < 0.05 from exploratory outcomes ##############
## Fractures PCA ########################################################################################
## Insert Status Fractures column in ds dataframe
ds <- data.frame(dataset['Status_Fractures'], norm_scale)
## Change Batch to categorical variable
ds$Status_Fractures <- factor(ds$Status_Fractures, levels=c(0,1), labels=c("No", "Yes"))

## PCA experiment - Set those variables which have shown a p-value < 0.05 with mann-whitney
# pca_dataset <- ds[ , c("Acetate", "Lactate", "Gly", "DHA", "GlycA", "Omega.3", "Unsaturation", "Ala", "Phosphoglyc",
#                        "Cholines","XL.HDL.L","M.HDL.FC","HDLsize","XL.HDL.PL","LA","L.HDL.PL","M.HDL.PL","HDL.PL")]
## PCA experiment - Set those variable with p-value < 0.05 and are on opposite sides of experiment above
pca_dataset <- ds[ , c("Gly", "Omega.3", "Ala","M.HDL.FC","L.HDL.PL","M.HDL.PL","HDL.PL")]

## PCA ##################################################################################################
pca <- prcomp(pca_dataset)
pca

## returns the scree plot showing the variances explained by the principal ## components versus the principal component number
plot(pca, type = "l")

# Get summary of the pca variable
summary(pca)  # Cutoff after 2/3

## Plot PCA #############################################################################################
## 2D PCA ###############################################################################################
library(ggfortify)
pca_dataset$Status_Fractures <- ds$Status_Fractures
autoplot(pca, data = pca_dataset, colour = 'Status_Fractures', label = FALSE)

library(ggbiplot)
g1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = ds$Status_Fractures, ellipse = TRUE, circle = TRUE)
g1 <- g1 + scale_color_discrete(name = '')
g1 <- g1 + theme(legend.direction = 'horizontal', legend.position = 'top')
g1
#########################################################################################################

## BMD PCA ##############################################################################################
## Insert Status BMB column in ds dataframe
ds <- data.frame(dataset['Status_BMD'], norm_scale)
## Change Status_BMD to categorical variable
ds$Status_BMD <- factor(ds$Status_BMD, levels=c(0,1,2), labels=c("Normal", "Osteopenic", "Osteoporotic"))

## PCA experiment - Get dataset to work on, i.e. remove Status_BMD which will be our variable
# pca_dataset <- ds[ , c("His", "Val", "Leu", "TotalBCAA", "Ile", "Phe", "Tyr")]
## PCA experiment with opposing effects on PCA
pca_dataset <- ds[ , c("His", "Gly", "Ala", "Cholines","Phosphoglyc")]

## PCA ##################################################################################################
pca <- prcomp(pca_dataset)
pca

## returns the scree plot showing the variances explained by the principal ## components versus the principal component number
plot(pca, type = "l")

# Get summary of the pca variable
summary(pca)  # Cutoff after 2

## Plot PCA #############################################################################################
## 2D PCA ###############################################################################################
library(ggfortify)
pca_dataset$Status_BMD <- ds$Status_BMD
autoplot(pca, data = pca_dataset, colour = 'Status_BMD', label = FALSE)

library(ggbiplot)
g1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = ds$Status_BMD, ellipse = TRUE, circle = TRUE)
g1 <- g1 + scale_color_discrete(name = '')
g1 <- g1 + theme(legend.direction = 'horizontal', legend.position = 'top')
g1
#########################################################################################################