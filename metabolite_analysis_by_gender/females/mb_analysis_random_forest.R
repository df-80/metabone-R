## IMPORTANT NOTE #######################################################################################
## Execute mb_data_filtering.R - From line 1 to 32 ######################################################
#########################################################################################################

## See: https://github.com/StatQuest/random_forest_demo/blob/master/random_forest_demo.R

## Load necessary libraries
pacman::p_load(cowplot, ggplot2, randomForest)

## Replace "TAG"s with NAs.
ds[ds == "TAG"] <- NA
# Remove columns which have more than 25 NAs
ds <- ds[, colSums(is.na(ds)) <= 24]

# Set fields to numeric
ds_metabolome <- sapply(ds, as.numeric)

## Adjust for Batch effect
ds_metabolome <- t(removeBatchEffect(t(ds_metabolome), batch = dataset$Batch))

###################################################################################################
## Status Fractures ###############################################################################
###################################################################################################

# Merge status fractures to the metabolome data
ds <- data.frame(dataset$Status_Fractures, ds_metabolome)
ds$dataset.Status_Fractures <- factor(ds$dataset.Status_Fractures, levels=c(0,1), labels = c("No", "Yes"))
names(ds)[names(ds) == "dataset.Status_Fractures"] <- "Status_Fractures"

# View dataset
View(ds)

# Show factors in the dataset
str(ds)

## Start building a random forest
set.seed(42)

## Impute any missing values in the training set using proximities
ds.imputed <- rfImpute(Status_Fractures ~ ., data = ds, iter = 6)

## Get random forest model with Status Fractures as the response variable
model <- randomForest(Status_Fractures ~ ., data=ds.imputed, proximity=TRUE)
## Print model
model

## Check if the random forest is big enough.
oob.error.data <- data.frame(
  Trees=rep(seq_len(nrow(model$err.rate)), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
    model$err.rate[,"Yes"],
    model$err.rate[,"No"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))

model <- randomForest(Status_Fractures ~ ., data=ds.imputed, ntree=1000, proximity=TRUE)
model

oob.error.data <- data.frame(
  Trees=rep(1:nrow(model$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
    model$err.rate[,"Yes"],
    model$err.rate[,"No"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))


oob.values <- vector(length=100)
for(i in 1:100) {
  temp.model <- randomForest(Status_Fractures ~ ., data=ds.imputed, mtry=i, ntree=1000)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate),1]
}
oob.values
## find the minimum error
min(oob.values)
## find the optimal value for mtry...
which(oob.values == min(oob.values))
## create a model for proximities using the best value for mtry
model <- randomForest(Status_Fractures ~ ., data=ds.imputed, ntree=1000, proximity=TRUE, mtry=which(oob.values == min(oob.values)))
model

## Start by converting the proximity matrix into a distance matrix.
distance.matrix <- as.dist(1-model$proximity)

mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=rownames(mds.values), X=mds.values[,1], Y=mds.values[,2], Status=ds.imputed$Status_Fractures)

ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(color=Status)) +
  theme_bw() +
  xlab(paste0("MDS1 - ", mds.var.per[1], "%")) +
  ylab(paste0("MDS2 - ", mds.var.per[2], "%")) +
  ggtitle("MDS plot using (1 - Random Forest Proximities)")


###################################################################################################
## Status BMD #####################################################################################
###################################################################################################

# Merge status BMD to the metabolome data
ds <- data.frame(dataset$Status_BMD, ds_metabolome)
ds$dataset.Status_BMD <- factor(ds$dataset.Status_BMD, levels=c(0,1,2), labels = c("Normal", "Osteopenia", "Osteoporosis"))
names(ds)[names(ds) == "dataset.Status_BMD"] <- "Status_BMD"

# View dataset
View(ds)

# Show factors in the dataset
str(ds)

## Start building a random forest
set.seed(42)

## Impute any missing values in the training set using proximities
ds.imputed <- rfImpute(Status_BMD ~ ., data = ds, iter=6)

## Get random forest model with Status Fractures as the response variable
model <- randomForest(Status_BMD ~ ., data=ds.imputed, proximity=TRUE)
## Print model
model

## Check if the random forest is big enough.
oob.error.data <- data.frame(
  Trees=rep(1:nrow(model$err.rate), times=4),
  Type=rep(c("OOB", "Normal", "Osteopenia", "Osteoporosis"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
    model$err.rate[,"Normal"],
    model$err.rate[,"Osteopenia"],
    model$err.rate[,"Osteoporosis"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))

model <- randomForest(Status_BMD ~ ., data=ds.imputed, ntree=200, proximity=TRUE)
model

oob.error.data <- data.frame(
  Trees=rep(1:nrow(model$err.rate), times=4),
  Type=rep(c("OOB", "Normal", "Osteopenia", "Osteoporosis"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
    model$err.rate[,"Normal"],
    model$err.rate[,"Osteopenia"],
    model$err.rate[,"Osteoporosis"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))


oob.values <- vector(length=10)
for(i in 1:10) {
  temp.model <- randomForest(Status_BMD ~ ., data=ds.imputed, mtry=i, ntree=1000)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate),1]
}
oob.values
## find the minimum error
min(oob.values)
## find the optimal value for mtry...
which(oob.values == min(oob.values))
## create a model for proximities using the best value for mtry
model <- randomForest(Status_BMD ~ ., data=ds.imputed, ntree=1000, proximity=TRUE, mtry=which(oob.values == min(oob.values)))
model

## Start by converting the proximity matrix into a distance matrix.
distance.matrix <- as.dist(1-model$proximity)

mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=rownames(mds.values), X=mds.values[,1], Y=mds.values[,2], Status=ds.imputed$Status_BMD)

ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(color=Status)) +
  theme_bw() +
  xlab(paste0("MDS1 - ", mds.var.per[1], "%")) +
  ylab(paste0("MDS2 - ", mds.var.per[2], "%")) +
  ggtitle("MDS plot using (1 - Random Forest Proximities)")