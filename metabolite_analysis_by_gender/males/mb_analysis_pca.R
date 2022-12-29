## IMPORTANT NOTE #######################################################################################
## Execute mb_data_filtering.R - From line 1 to 63 ######################################################
## To see PCA before batch effect correction, execute from line 1 to line 45 ############################
#########################################################################################################

# Load necessary libraries
pacman::p_load(caret, ggplot2)

# Filter by Gender - get Males only
ds <- ds[ds$Gender == 'Male',]
cat("Rows:", nrow(ds),", Columns:", ncol(ds),"\n")

# Save the non-metabolite variables
nmv_ds <- ds[ , 1:22]
ds <- ds[ , 23:ncol(ds)]

## Scale data - Use caret library to MinMax scale all columns
process <- preProcess(as.data.frame(ds), method = "range")
norm_scale <- predict(process, as.data.frame(ds))


## Fractures PCA ########################################################################################
## Insert Status Fractures column in ds dataframe
ds <- data.frame(nmv_ds['Status_Fractures'], norm_scale)

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
ds <- data.frame(nmv_ds['Status_BMD'], norm_scale)

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
ds <- data.frame(nmv_ds['LS_status'], norm_scale)

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
ds <- data.frame(nmv_ds['FN_status'], norm_scale)

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
ds <- data.frame(nmv_ds['TH_status'], norm_scale)

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
ds <- data.frame(nmv_ds['Status_Fractures'], norm_scale)

## PCA experiment - Set those variables which have shown a p-value < 0.05 with mann-whitney
pca_dataset <- ds[ , c("Gly", "Total.CE", "Sphingomyelins", "Total.C", "Total.FC", "ApoB", "Total.L", "Omega.6")]
# pca_dataset <- ds[ , c("Gly", "L.LDL.FC", "LDL.FC", "M.LDL.FC", "Clinical.LDL.C", "Total.CE", "L.LDL.PL", "L.LDL.C", "L.LDL.CE", "S.LDL.FC",
#                        "Sphingomyelins", "LDL.PL", "L.LDL.L", "LDL.C", "S.LDL.C", "M.VLDL.CE", "S.LDL.P", "Total.C", "LDL.CE", "LDL.L",
#                        "IDL.C", "IDL.FC", "L.LDL.P", "M.LDL.C", "S.LDL.L", "S.LDL.CE", "Total.FC", "S.LDL.PL", "IDL.CE", "M.LDL.PL",
#                        "non.HDL.C", "M.VLDL.C", "M.LDL.L", "M.LDL.CE", "IDL.L", "IDL.PL", "LDL.P", "IDL.P", "ApoB", "Remnant.C", "M.HDL.TG",
#                        "M.VLDL.FC", "S.VLDL.FC", "Total.L", "Omega.6", "HDL.TG", "M.LDL.P")]
# pca_dataset <- ds[ , c("Gly", "LDL.FC", "Clinical.LDL.C", "Total.CE", "Sphingomyelins", "LDL.PL", "LDL.C", "Total.C", "LDL.CE", "LDL.L",
#                        "IDL.C", "IDL.FC", "Total.FC", "IDL.CE", "non.HDL.C", "IDL.L", "IDL.PL", "LDL.P", "IDL.P", "ApoB", "Remnant.C",
#                        "Total.L", "Omega.6", "HDL.TG")]

## PCA ##################################################################################################
pca <- prcomp(pca_dataset)
pca

## returns the scree plot showing the variances explained by the principal components
## versus the principal component number
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
## Insert Status BMD column in ds dataframe
ds <- data.frame(nmv_ds['Status_BMD'], norm_scale)

## PCA experiment - Set those variables which have shown a p-value < 0.05 with kruskall-wallis
pca_dataset <- ds[ , c("Sphingomyelins", "Gly", "Total.C", "ApoB", "Total.CE", "His", "Total.L")]

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

## Total Hip BMD PCA ####################################################################################
## Insert Status BMD column in ds dataframe
ds <- data.frame(nmv_ds['TH_status'], norm_scale)

## PCA experiment - Set those variables which have shown a p-value < 0.05 with kruskall-wallis
pca_dataset <- ds[ , c("Gly", "Acetoacetate", "Acetone", "Sphingomyelins","LDL.PL")]

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
pca_dataset$TH_status <- ds$TH_status
autoplot(pca, data = pca_dataset, colour = 'TH_status', label = FALSE)

library(ggbiplot)
g1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = ds$TH_status, ellipse = TRUE, circle = TRUE)
g1 <- g1 + scale_color_discrete(name = '')
g1 <- g1 + theme(legend.direction = 'horizontal', legend.position = 'top')
g1
#########################################################################################################