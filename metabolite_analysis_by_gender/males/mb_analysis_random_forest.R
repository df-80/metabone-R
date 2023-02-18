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
ds_metabolome <- data.frame(dataset$Gender, dataset$Status_Fractures, dataset$Age, dataset$Alcohol, dataset$BMI, ds)
ds_metabolome$dataset.Status_Fractures <- factor(ds_metabolome$dataset.Status_Fractures, levels=c(0,1), labels = c("No", "Yes"))
names(ds_metabolome)[names(ds_metabolome) == "dataset.Gender"] <- "Gender"
names(ds_metabolome)[names(ds_metabolome) == "dataset.Status_Fractures"] <- "Status_Fractures"
names(ds_metabolome)[names(ds_metabolome) == "dataset.Age"] <- "Age"
names(ds_metabolome)[names(ds_metabolome) == "dataset.Alcohol"] <- "Alcohol"
names(ds_metabolome)[names(ds_metabolome) == "dataset.BMI"] <- "BMI"

# Filter by Gender - get Males only
ds_metabolome <- ds_metabolome[ds_metabolome$Gender == 0, ]

# Drop Gender columns as all values are 0.
ds_metabolome <- ds_metabolome[-1]

# View dataset
View(ds_metabolome)

# Show factors in the dataset
str(ds_metabolome)

## Start building a random forest
set.seed(42)

## Impute any missing values in the training set using proximities - No need to impute for Males
#ds.imputed <- rfImpute(Status_Fractures ~ ., data = ds_metabolome, iter = 6)

## Get random forest model with Status Fractures as the response variable
#model <- randomForest(Status_Fractures ~ ., data=ds_metabolome, proximity=TRUE)
model <- randomForest(
  Status_Fractures ~ Gly + L.LDL.FC + LDL.FC + M.LDL.FC + Clinical.LDL.C + Total.CE + L.LDL.PL + L.LDL.C + L.LDL.CE + S.LDL.FC +
    Sphingomyelins + LDL.PL + L.LDL.L + LDL.C + S.LDL.C + M.VLDL.CE + S.LDL.P + Total.C + LDL.CE + LDL.L + IDL.C + IDL.FC + L.LDL.P +
    M.LDL.C + S.LDL.L + S.LDL.CE + Total.FC + S.LDL.PL + IDL.CE + M.LDL.PL + non.HDL.C + M.VLDL.C + M.LDL.L + M.LDL.CE + IDL.L + IDL.PL +
    LDL.P + IDL.P + ApoB + Remnant.C + M.HDL.TG + M.VLDL.FC + S.VLDL.FC + Total.L + Omega.6 + HDL.TG + M.LDL.P,
  data=ds_metabolome, proximity=TRUE)
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

# model <- randomForest(Status_Fractures ~ ., data=ds.imputed, ntree=1000, proximity=TRUE)
# model
#
# oob.error.data <- data.frame(
#   Trees=rep(1:nrow(model$err.rate), times=3),
#   Type=rep(c("OOB", "Yes", "No"), each=nrow(model$err.rate)),
#   Error=c(model$err.rate[,"OOB"],
#     model$err.rate[,"Yes"],
#     model$err.rate[,"No"]))
#
# ggplot(data=oob.error.data, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))


oob.values <- vector(length=100)
for(i in 1:100) {
  temp.model <- randomForest(Status_Fractures ~ ., data=ds_metabolome, mtry=i, ntree=1000)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate),1]
}
oob.values
## find the minimum error
min(oob.values)
## find the optimal value for mtry...
which(oob.values == min(oob.values))
## create a model for proximities using the best value for mtry
model <- randomForest(Status_Fractures ~ ., data=ds_metabolome, ntree=1000, proximity=TRUE, mtry=which(oob.values == min(oob.values)))
model

## Start by converting the proximity matrix into a distance matrix.
distance.matrix <- as.dist(1-model$proximity)

mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=rownames(mds.values), X=mds.values[,1], Y=mds.values[,2], Status=ds_metabolome$Status_Fractures)

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

## Impute any missing values in the training set using proximities - No need to impute for Males
# ds.imputed <- rfImpute(Status_BMD ~ ., data = ds, iter=6)

## Get random forest model with Status Fractures as the response variable
#model <- randomForest(Status_BMD ~ ., data=ds.imputed, proximity=TRUE)
model <- randomForest(
  Status_BMD ~ Sphingomyelins + Gly + M.VLDL.FC + S.LDL.P + Total.C + S.LDL.CE + LDL.P + M.VLDL.CE + M.VLDL.C + L.LDL.P + ApoB + S.LDL.L +
    S.LDL.PL + Glycerol + Glucose + S.LDL.C + non.HDL.C + L.LDL.CE + IDL.P + S.VLDL.FC + L.LDL.C + M.VLDL.PL + L.LDL.FC + Total.CE +
    Clinical.LDL.C + HDLsize + M.LDL.P + His + Total.L + L.LDL.L + IDL.FC + LDL.CE + LDL.L + LDL.C + M.LDL.CE + LDL.PL + L.LDL.PL +
    M.LDL.L + M.VLDL.P,
  data=ds_metabolome, proximity=TRUE)

## Print model
model

## Check if the random forest is big enough.
oob.error.data <- data.frame(
  Trees=rep(seq_len(nrow(model$err.rate)), times=4),
  Type=rep(c("OOB", "Normal", "Osteopenia", "Osteoporosis"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
    model$err.rate[,"Normal"],
    model$err.rate[,"Osteopenia"],
    model$err.rate[,"Osteoporosis"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))

model <- randomForest(Status_BMD ~ ., data=ds_metabolome, ntree=200, proximity=TRUE)
model

oob.error.data <- data.frame(
  Trees=rep(seq_len(nrow(model$err.rate)), times=4),
  Type=rep(c("OOB", "Normal", "Osteopenia", "Osteoporosis"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
    model$err.rate[,"Normal"],
    model$err.rate[,"Osteopenia"],
    model$err.rate[,"Osteoporosis"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))


oob.values <- vector(length=10)
for(i in 1:10) {
  temp.model <- randomForest(Status_BMD ~ ., data=ds_metabolome, mtry=i, ntree=1000)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate),1]
}
oob.values
## find the minimum error
min(oob.values)
## find the optimal value for mtry...
which(oob.values == min(oob.values))
## create a model for proximities using the best value for mtry
model <- randomForest(Status_BMD ~ ., data=ds_metabolome, ntree=1000, proximity=TRUE, mtry=which(oob.values == min(oob.values)))
model

## Start by converting the proximity matrix into a distance matrix.
distance.matrix <- as.dist(1-model$proximity)

mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=rownames(mds.values), X=mds.values[,1], Y=mds.values[,2], Status=ds_metabolome$Status_BMD)

ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
  geom_text(aes(color=Status)) +
  theme_bw() +
  xlab(paste0("MDS1 - ", mds.var.per[1], "%")) +
  ylab(paste0("MDS2 - ", mds.var.per[2], "%")) +
  ggtitle("MDS plot using (1 - Random Forest Proximities)")