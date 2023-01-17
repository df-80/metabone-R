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
  norm_scale['Glucose'], norm_scale['Acetate'], norm_scale['XL.HDL.L'], norm_scale['Lactate'], norm_scale['GlycA'],
  norm_scale['XL.HDL.PL'], norm_scale['XL.HDL.C'], norm_scale['XL.HDL.FC'], norm_scale['XL.HDL.CE'],
  norm_scale['His'], norm_scale['XL.HDL.P'], norm_scale['Gly'])

# Filter by Gender - get Males only
ds_new <- ds_new[ds_new$Gender == 1, ]

# Drop Gender column as all values are 0.
ds_new <- ds_new[-2]

# Change Status Fractures to categorical variable
ds_new$Status_Fractures <- factor(ds_new$Status_Fractures, levels=c(0,1), labels=c("No", "Yes"))
cat("Rows:", nrow(ds_new),", Columns:", ncol(ds_new),"\n")
summary(ds_new$Status_Fractures)

# Remove columns which have more than 25 NAs
ds_new <- ds_new[, colSums(is.na(ds_new)) <= 24]
ds.noNA <- ds_new %>% na.omit()

cat("Rows:", nrow(ds.noNA),", Columns:", ncol(ds.noNA),"\n")
summary(ds.noNA$Status_Fractures)

balanced.ds <- ds.noNA %>% group_by(Status_Fractures) %>% slice_sample(n = 87)

## 70% of the sample size
train_size <- floor(0.7 * nrow(balanced.ds))

set.seed(55)
in_rows <- sample(seq_len(nrow(balanced.ds)), size = train_size, replace = FALSE)

train <- balanced.ds[in_rows, ]
test  <- balanced.ds[-in_rows, ]
attach(train)

mult_result <- glm(Status_Fractures ~
                     Age + BMI + Acetate + XL.HDL.L + Lactate + GlycA + XL.HDL.PL + XL.HDL.C + XL.HDL.FC + XL.HDL.CE + His + XL.HDL.P + Gly,
                   family = 'binomial', data=train)
mult_result

predict <- predict(mult_result, newdata=test, type="response")
predict.cat <- cut(predict,breaks=c(0, 0.5, 1), labels=c("No fractures", "Fractures"))
summary(predict.cat)
summary(test$Status_Fractures)
table(predict.cat, test$Status_Fractures)

svmfit  <- svm(Status_Fractures ~
#                 Age + BMI +
                 Acetate + XL.HDL.L + Lactate + GlycA + XL.HDL.PL + XL.HDL.C + XL.HDL.FC + XL.HDL.CE + His + XL.HDL.P,
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
  norm_scale['XXL.VLDL.TG'], norm_scale['Val'], norm_scale['Citrate'], norm_scale['Leu'], norm_scale['TotalBCAA'],
  norm_scale['His'], norm_scale['XXL.VLDL.PL'], norm_scale['XXL.VLDL.L'], norm_scale['VLDLsize'], norm_scale['XXL.VLDL.P'],
  norm_scale['Ile'], norm_scale['XXL.VLDL.C'], norm_scale['XXL.VLDL.FC'], norm_scale['Phe'])

# Filter by Gender - get Males only
ds_new <- ds_new[ds_new$Gender == 1, ]

# Drop Gender column as all values are 0.
ds_new <- ds_new[-2]

# Change BMD to categorical variable
ds_new$Status_BMD <- factor(ds_new$Status_BMD, levels=c(0,1,2), labels=c("Normal", "Osteopenic", "Osteoporotic"))
cat("Rows:", nrow(ds_new),", Columns:", ncol(ds_new),"\n")
summary(ds_new$Status_BMD)

# Remove columns which have more than 25 NAs
ds_new <- ds_new[, colSums(is.na(ds_new)) <= 24]
ds.noNA <- ds_new %>% na.omit()

cat("Rows:", nrow(ds.noNA),", Columns:", ncol(ds.noNA),"\n")
summary(ds.noNA$Status_BMD)

balanced.ds <- ds.noNA %>% group_by(Status_BMD) %>% slice_sample(n = 36)

## 70% of the sample size
train_size <- floor(0.7 * nrow(balanced.ds))

set.seed(33)
in_rows <- sample(seq_len(nrow(balanced.ds)), size = train_size, replace = FALSE)

train <- balanced.ds[in_rows, ]
test  <- balanced.ds[-in_rows, ]
attach(train)

## Multinomial Regression #############################################################################
mult_result <- multinom(train$Status_BMD ~
                          Age + BMI + XXL.VLDL.TG + Val + Leu + TotalBCAA + His + XXL.VLDL.PL +
                            XXL.VLDL.L + VLDLsize + XXL.VLDL.P + Ile + XXL.VLDL.C + XXL.VLDL.FC + Phe,
                        maxit = 100000, data=train)

summary(mult_result)
tab_model(mult_result, show.est=FALSE)

predict <- predict(mult_result, newdata=test, type="class")
actual <- test$Status_BMD
confusionMatrix(table(actual, predict))

## SVM ################################################################################################
svmfit <- svm(train$Status_BMD ~
                Age + BMI + Val + Leu + TotalBCAA + His + XXL.VLDL.PL + XXL.VLDL.L + VLDLsize + Ile + Phe,
               data = train, kernel='linear', cost = 10, scale = TRUE)

predict <- predict(svmfit, newdata=test, type="response")
actual  <- test$Status_BMD
confusionMatrix(table(actual, predict))