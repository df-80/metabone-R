## IMPORTANT NOTE #######################################################################################
## Execute mb_data_filtering.R - From line 1 to 32 and 44 to 51 #########################################
#########################################################################################################

# Load necessary libraries
pacman::p_load(caret, e1071, nnet, ggplot2, sjPlot)

ds <- data.frame(dataset['Age'], dataset['BMI'], ds)

## Scale data - Use caret library to MinMax scale all columns
process <- preProcess(as.data.frame(ds), method = "range")
norm_scale <- predict(process, as.data.frame(ds))

###################################################################################################
## Status Fractures ###############################################################################
###################################################################################################

# Add Status_Fractures, Gender, Smoking, Alcohol, Age and BMI to a new dataframe
ds_new <- data.frame(
  dataset['Status_Fractures'], dataset['Gender'], dataset['Smoking'],  dataset['Alcohol'], norm_scale['Age'], norm_scale['BMI'],
  norm_scale['Acetate'], norm_scale['Lactate'], norm_scale['Glucose'], norm_scale['Gly'], norm_scale['DHA'],
  norm_scale['GlycA'], norm_scale['Phosphatidylc'], norm_scale['Omega.3'], norm_scale['Unsaturation'],
  norm_scale['Ala'], norm_scale['Phosphoglyc'], norm_scale['Cholines'], norm_scale['XL.HDL.L'], norm_scale['Glycerol'],
  norm_scale['M.HDL.FC'], norm_scale['HDLsize'], norm_scale['XL.HDL.PL'], norm_scale['LA'])

# Change Status Fractures to categorical variable
ds_new$Status_Fractures <- factor(ds_new$Status_Fractures, levels=c(0,1), labels=c("No", "Yes"))

# Remove columns which have more than 25 NAs
ds_new <- ds_new[, colSums(is.na(ds_new)) <= 24]
ds.noNA <- ds_new %>% na.omit()

summary(ds.noNA)
balanced.ds <- ds.noNA %>% group_by(Status_Fractures) %>% slice_sample(n = 92)

## 70% of the sample size
train_size <- floor(0.7 * nrow(balanced.ds))

set.seed(55)
in_rows <- sample(seq_len(nrow(balanced.ds)), size = train_size, replace = FALSE)

train <- balanced.ds[in_rows, ]
test  <- balanced.ds[-in_rows, ]
attach(train)

mult_result <- glm(Status_Fractures ~ Age + BMI + Acetate + Lactate + Gly + DHA + GlycA + Phosphatidylc + Omega.3 +
                       Unsaturation + Ala + Phosphoglyc + Cholines + XL.HDL.L + M.HDL.FC + HDLsize + XL.HDL.PL + LA,
                   family = 'binomial', data=train)
mult_result

predict <- predict(mult_result, newdata=test, type="response")
predict.cat <- cut(predict,breaks=c(0, 0.5, 1), labels=c("No fractures", "Fractures"))
summary(predict.cat)
summary(test$Status_Fractures)
table(predict.cat, test$Status_Fractures)

svmfit  <- svm(train$Status_Fractures ~ Age + BMI + Acetate + Lactate + Gly + DHA + GlycA + Phosphatidylc + Omega.3 +
                       Unsaturation + Ala + Phosphoglyc + Cholines + XL.HDL.L + M.HDL.FC + HDLsize + XL.HDL.PL + LA,
               data = train, kernel = "linear", cost = 10, scale = FALSE)

predict <- predict(svmfit, newdata=test, type="response")
actual  <- test$Status_Fractures
confusionMatrix(table(predict, actual))


###################################################################################################
## Status BMD #####################################################################################
###################################################################################################

# Add Status_Fractures, Gender, Smoking, Alcohol, Age and BMI to a new dataframe
ds_new <- data.frame(
  dataset['Status_BMD'], dataset['Gender'], dataset['Smoking'],  dataset['Alcohol'], norm_scale['Age'], norm_scale['BMI'],
  norm_scale['His'], norm_scale['Val'], norm_scale['Citrate'], norm_scale['Leu'], norm_scale['TotalBCAA'], norm_scale['Ile'], norm_scale['Phe'])

# Change BMD to categorical variable
ds_new$Status_BMD <- factor(ds_new$Status_BMD, levels=c(0,1,2), labels=c("Normal", "Osteopenic", "Osteoporotic"))

# Remove columns which have more than 25 NAs
ds.noNA <- ds_new %>% na.omit()

summary(ds.noNA)
balanced.ds <- ds.noNA %>% group_by(Status_BMD) %>% slice_sample(n = 22)

## 70% of the sample size
train_size <- floor(0.7 * nrow(balanced.ds))

set.seed(33)
in_rows <- sample(seq_len(nrow(balanced.ds)), size = train_size, replace = FALSE)

train <- balanced.ds[in_rows, ]
test  <- balanced.ds[-in_rows, ]
attach(train)

## Multinomial Regression #############################################################################
mult_result <- multinom(
  train$Status_BMD ~ His,
  maxit = 100000, data=train)

summary(mult_result)
tab_model(mult_result, show.est=FALSE)

predict <- predict(mult_result, newdata=test, type="class")
actual <- test$Status_BMD
confusionMatrix(table(predict, actual))

## SVM ################################################################################################
svmfit  <- svm(
  train$Status_BMD ~
    His + Citrate, data = train, kernel='radial', cost = 10, scale = TRUE)
predict <- predict(svmfit, newdata=test, type="response")
actual  <- test$Status_BMD
confusionMatrix(table(predict, actual))

###################################################################################################
## Lumbar Spine ###################################################################################
###################################################################################################

# Add Status_Fractures, Gender, Smoking, Alcohol, Age and BMI to a new dataframe
ds_new <- data.frame(
  dataset['LS_status'], dataset['Gender'], dataset['Smoking'],  dataset['Alcohol'], norm_scale['Age'], norm_scale['BMI'],
  norm_scale['LDLsize'], norm_scale['Citrate'], norm_scale['His'])

# Change Lumbar Spine status to categorical variable
ds_new$LS_status <- factor(ds_new$LS_status, levels=c(0,1,2), labels=c("Normal", "Osteopenic", "Osteoporotic"))

# Remove columns which have more than 25 NAs
ds.noNA <- ds_new %>% na.omit()

summary(ds.noNA)
balanced.ds <- ds.noNA %>% group_by(LS_status) %>% slice_sample(n = 44)

## 70% of the sample size
train_size <- floor(0.7 * nrow(balanced.ds))

set.seed(47)
in_rows <- sample(seq_len(nrow(balanced.ds)), size = train_size, replace = FALSE)

train <- balanced.ds[in_rows, ]
test  <- balanced.ds[-in_rows, ]
attach(train)

## Multinomial Regression #############################################################################
mult_result <- multinom(
  train$LS_status ~ Citrate + His,
  maxit = 100000, data=train)

summary(mult_result)
tab_model(mult_result, show.est=FALSE)

predict <- predict(mult_result, newdata=test, type="class")
actual <- test$LS_status
confusionMatrix(table(predict, actual))

## SVM ################################################################################################
svmfit  <- svm(
  train$LS_status ~
    Citrate + His, data = train, kernel='radial', cost = 10, scale = TRUE)
predict <- predict(svmfit, newdata=test, type="response")
actual  <- test$LS_status
confusionMatrix(table(predict, actual))