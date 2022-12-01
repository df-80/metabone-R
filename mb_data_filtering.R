## Remove all variables from the environment
rm(list = ls())

## Install pacman package if not available. Pacman is a package
## management tool which installs and loads the necessary libraries.
if (!require("pacman"))
  install.packages("pacman")

## Load necessary libraries
pacman::p_load(pacman, limma, dplyr, DMwR2)

## Import dataset #######################################################################################
dataset <- read.csv(file.choose(), header = TRUE)
cat("Rows:", nrow(dataset),", Columns:", ncol(dataset),"\n")

## Define functions #####################################################################################
# define function to count number of NAs
cnt_na <- apply(dataset, 1, function(z) sum(is.na(z)))


## Data Filtering #######################################################################################
# Remove rows with many NAs (>= 50 columns with NAs)
dataset <- dataset[cnt_na < 50, ]
cat("Rows:", nrow(dataset),", Columns:", ncol(dataset),"\n")

# Get from 23rd to 193st column from dataset dataframe
ds <- as.data.frame(dataset[, 23:193])
# Remove TG.PG and ApoB.ApoA1 columns
ds <- ds[, -c(which(colnames(ds)=="TG.PG"), which(colnames(ds)=="ApoB.ApoA1"))]
# Count number of rows and columns
cat("Rows:", nrow(ds),", Columns:", ncol(ds),"\n")

## Optional steps for certain data analysis algorithms
# Remove columns which have more than 25 NAs
ds <- ds[, colSums(is.na(ds)) <= 24]
# Remove columns with Tags
ds %>% select_if(~any(grepl("TAG", .)))
ds <- ds[, -c(which(colnames(ds) == "Glycerol"), which(colnames(ds) == "bOHbutyrate"))]
# Impute NAs
ds <- knnImputation(ds, k = 10, scale = FALSE, meth = "weighAvg", distData = NULL)

# Set all columns to numeric
ds <- sapply(ds, as.numeric)
# Count number of rows and columns and view data
cat("Rows:", nrow(ds),", Columns:", ncol(ds),"\n")
View(ds)

## Adjust for Batch effect
# Transpose the dataframe and return a matrix - t(ds)
ds <- as.data.frame(t(removeBatchEffect(t(ds), batch = dataset$Batch)))

## Merge first 22 columns of original dataset to ds
ds <- cbind(dataset[1:22], ds)

## Set categorical variables
ds$FS <- factor(ds$FS, levels = c(0,1), labels = c("No", "Yes"))
ds$Gender <- factor(ds$Gender,  levels=c(0,1), labels=c("Male", "Female"))
ds$Status_BMD <- factor(ds$Status_BMD, levels=c(0,1,2), labels=c("Normal", "Osteopenia", "Osteoporosis"))
ds$LS_status <- factor(ds$LS_status, levels=c(0,1,2), labels=c("Normal", "Osteopenia", "Osteoporosis"))
ds$FN_status <- factor(ds$FN_status, levels=c(0,1,2), labels=c("Normal", "Osteopenia", "Osteoporosis"))
ds$TH_status <- factor(ds$TH_status, levels=c(0,1,2), labels=c("Normal", "Osteopenia", "Osteoporosis"))
ds$Status_Fractures <- factor(ds$Status_Fractures, levels=c(0,1), labels = c("No", "Yes"))

## Attach table to access columns directly by name
attach(ds)

## Summary of a number of variables
summary(ds$Gender)
summary(ds$Age)
summary(ds$Status_Fractures)
summary(ds$Status_BMD)
summary(ds$LS_status)
summary(ds$FN_status)
summary(ds$TH_status)