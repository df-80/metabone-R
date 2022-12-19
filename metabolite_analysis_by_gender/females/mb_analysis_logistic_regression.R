## IMPORTANT NOTE #######################################################################################
## Execute mb_data_filtering.R - From line 1 to 57 and 59 to 63 #########################################
#########################################################################################################

# Load necessary libraries
pacman::p_load(caret, nnet, sjPlot)

# Setting explanatory variables only
additional_columns <- ds[, c("Gender", "Age", "BMI", "Smoking", "Alcohol")]
ds_new <- cbind(additional_columns, ds[, 23:ncol(ds)])

head(ds_new)
attach(ds_new)

### Find which attributes are highly correlated - multicolliniarity analysis
# calculate correlation matrix
correlationMatrix <- round(cor(ds_new, method='spearman'), 2)
# find attributes that are highly corrected
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
print(highlyCorrelated)
# Not hightly correlated @ 0.5 - 1 - 6, 38, 47, 57, 65-68, 141

######## Calculate Status BMD from manually selecting the variable with lower correlation
ds_new_status_bmd <- cbind(ds_new, Status_BMD)

mult_result <- multinom(ds_new_status_bmd$Status_BMD ~
                          Gender + Age + BMI + Smoking + Alcohol + Total.C + LDLsize + Unsaturation + His + Acetate + Acetoacetate + Acetone + Albumin + XL.HDL.L,
                        maxit = 1000, data = ds_new_status_bmd)
summary(mult_result)

## Reduced the variables to 7 according to the standard error vs the model coefficients - STD err has to be less than half of model co-off
mult_result <- multinom(ds_new_status_bmd$Status_BMD ~ Age + BMI + Smoking + LDLsize + His + Acetoacetate + Acetone,
                        maxit = 1000, data = ds_new_status_bmd)
summary(mult_result)
tab_model(mult_result, show.est=FALSE)

# Test whether the variables have an effect by testing the intercept only
multi_intercept_only <- multinom(ds_new_status_bmd$Status_BMD ~ 1, data=ds_new_status_bmd)
summary(multi_intercept_only)

# With a p-value (Pr(Chi)) of 0 < 0.05, the change in deviance is significant.
# Hence there is an improvement in the model with the variables.
anova(multi_intercept_only, mult_result, test='Chisq')

# Prediction estimates
predicted<-predict(mult_result, ds_new_status_bmd$Status_BMD, type = "probs")
head(predicted)
predicted<-predict(mult_result, ds_new_status_bmd$Status_BMD, type = "class")
head(predicted)

# Convert prediction to classifications
conttable <- with(ds_new_status_bmd,table(predicted, Status_BMD))
conttable

# Calculate the percentages of the predictions
percentCorrect <- round(diag(conttable)/rowSums(conttable),2)
finaltable <- cbind(conttable, percentCorrect)
finaltable
#########################################################################################################

######## Calculate Fracture Status from manually selecting the variable with lower correlation
ds_new_status_fractures <- cbind(ds_new, Status_Fractures)

lm_result <- glm(ds_new_status_fractures$Status_Fractures ~
                          Gender + Age + BMI + Smoking + Alcohol + Total.C + LDLsize + Unsaturation + His + Acetate + Acetoacetate + Acetone + Albumin + XL.HDL.L,
                        family = 'binomial', data = ds_new_status_fractures)
summary(lm_result)

## Reduced the variables to 5 according to the standard error vs the model coefficients - STD err has to be less than half of model co-off
lm_result <- glm(ds_new_status_fractures$Status_Fractures ~ Gender + His + Acetate + Acetone + Gly,
                        family = 'binomial', data = ds_new_status_fractures)
summary(lm_result)
tab_model(lm_result, show.est=FALSE)

# Test whether the variables have an effect by testing the intercept only
multi_intercept_only <- glm(ds_new_status_fractures$Status_Fractures ~ 1, family = 'binomial', data=ds_new_status_fractures)
summary(multi_intercept_only)

# With a p-value (Pr(Chi)) of 0 < 0.05, the change in deviance is significant.
# Hence there is an improvement in the model with the variables.
anova(multi_intercept_only, lm_result, test='Chisq')

# Prediction estimates
predicted <- plogis(predict(lm_result, ds_new_status_fractures))
head(predicted)

Pearson <- residuals(lm_result, type='pearson')
plot(lm_result$linear.predictors, Pearson)

which(abs(Pearson) > 3)
length(which(predicted == 0.5))

predicted.cat <- cut(predicted,breaks=c(0, 0.5, 1), labels=c("No fractures", "Fractures"))
summary(predicted.cat)
summary(ds_new_status_fractures$Status_Fractures)
table(predicted.cat, ds_new_status_fractures$Status_Fractures)
#########################################################################################################


### Use Caret library to automatically determine the features to select #################################
### BMD Feature Selection ###############################################################################
### Reference - https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
# Ensure results are repeatable
set.seed(42)

# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10)
# run the RFE algorithm - Status BMD
results_status_bmd <- rfe(ds_new, ds$Status_BMD, sizes = seq_len(ncol(ds_new)), rfeControl=control)
# summarize top 5 results - BMI, DHA, Omega.3, Age, S.HDL.TG
print(results_status_bmd)
# list the chosen features
predictors(results_status_bmd)
# plot the results
plot(results_status_bmd, type=c("g", "o"))

mult_result <- multinom(ds$Status_BMD ~ BMI + Age + Leu, maxit = 10000, data=ds_new)

summary(mult_result)
tab_model(mult_result, show.est=FALSE)

# Test whether the variables have an effect by testing the intercept only
multi_intercept_only <- multinom(ds$Status_BMD ~ 1, data=ds_new)
summary(multi_intercept_only)

# With a p-value (Pr(Chi)) of 0 < 0.05, the change in deviance is significant.
# Hence there is an improvement in the model with the variables.
anova(multi_intercept_only, mult_result, test='Chisq')

# Prediction estimates
predicted<-predict(mult_result, ds$Status_BMD, type = "probs")
head(predicted)
predicted<-predict(mult_result, ds$Status_BMD, type = "class")
head(predicted)

# Convert prediction to classifications
conttable <- with(ds_new,table(predicted, ds$Status_BMD))
conttable

# Calculate the percentages of the predictions
percentCorrect <- round(diag(conttable)/rowSums(conttable),2)
finaltable <- cbind(conttable, percentCorrect)
finaltable
#########################################################################################################


### Fractures Feature Selection #########################################################################
# run the RFE algorithm - Status Fractures
set.seed(47)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
results_status_fractures <- rfe(ds_new, ds$Status_Fractures, sizes = seq_len(ncol(ds_new)), rfeControl=control)
# summarize tp 5 results - Acetate, Gly, GlycA, Age, DHA
print(results_status_fractures)
# list the chosen features
predictors(results_status_fractures)
# plot the results
plot(results_status_fractures, type=c("g", "o"))

glr_result <- glm(ds$Status_Fractures ~ Acetate + GlycA + Age + DHA + XL.HDL.FC + Phosphatidylc + Ile,
                  data=ds_new, family = "binomial")

summary(glr_result)
tab_model(glr_result, show.est=FALSE)

# Test whether the variables have an effect by testing the intercept only
glr_intercept_only <- glm(ds$Status_Fractures ~ 1, data=ds_new, family = "binomial")
summary(glr_intercept_only)

# With a p-value (Pr(Chi)) of 0 < 0.05, the change in deviance is significant.
# Hence there is an improvement in the model with the variables.
anova(glr_intercept_only, glr_result, test='Chisq')

# Prediction estimates
predicted <- plogis(predict(glr_result, ds_new))
head(predicted)

Pearson <- residuals(glr_result, type='pearson')
plot(glr_result$linear.predictors, Pearson)

which(abs(Pearson) > 3)
length(which(predicted == 0.5))

predicted.cat <- cut(predicted,breaks=c(0, 0.5, 1), labels=c("No fractures", "Fractures"))
summary(predicted.cat)
summary(ds$Status_Fractures)
table(predicted.cat, ds$Status_Fractures)
#########################################################################################################


### Lumbar Spine Feature Selection ######################################################################
# run the RFE algorithm - Lumbar Spine
set.seed(43)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
results_lumbar_spine <- rfe(ds_new, ds$LS_status, sizes = seq_len(ncol(ds_new)), rfeControl=control)
# summarize top 5 results - BMI, XL.HDL.L, IDL.TG, Unsaturation, Gender
print(results_lumbar_spine)
# list the chosen features
predictors(results_lumbar_spine)
# plot the results
plot(results_lumbar_spine, type=c("g", "o"))

mult_result <- multinom(ds$LS_status ~ BMI + IDL.TG, maxit = 10000, data=ds_new)

summary(mult_result)
tab_model(mult_result, show.est=FALSE)

# Test whether the variables have an effect by testing the intercept only
multi_intercept_only <- multinom(ds$LS_status ~ 1, data=ds_new)
summary(multi_intercept_only)

# With a p-value (Pr(Chi)) of 0 < 0.05, the change in deviance is significant.
# Hence there is an improvement in the model with the variables.
anova(multi_intercept_only, mult_result, test='Chisq')

# Prediction estimates
predicted<-predict(mult_result, ds$LS_status, type = "probs")
head(predicted)
predicted<-predict(mult_result, ds$LS_status, type = "class")
head(predicted)

# Convert prediction to classifications
conttable <- with(ds_new, table(predicted, ds$LS_status))
conttable

# Calculate the percentages of the predictions
percentCorrect <- round(diag(conttable)/rowSums(conttable),2)
finaltable <- cbind(conttable, percentCorrect)
finaltable
#########################################################################################################


### Femoral Neck Feature Selection ######################################################################
# run the RFE algorithm - Femoral Neck
set.seed(7)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
results_femoral_neck <- rfe(ds_new, ds$FN_status, sizes = seq_len(ncol(ds_new)), rfeControl=control)
# summarize the results - Age, GlycA, XL.HDL.C, BMI, Acetoacetate
print(results_femoral_neck)
# list the chosen features
predictors(results_femoral_neck)
# plot the results
plot(results_femoral_neck, type=c("g", "o"))

mult_result <- multinom(ds$FN_status ~ Age + GlycA + XL.HDL.C + BMI + Acetoacetate, maxit = 10000, data=ds_new)

summary(mult_result)
tab_model(mult_result, show.est=FALSE)

# Test whether the variables have an effect by testing the intercept only
multi_intercept_only <- multinom(ds$FN_status ~ 1, data=ds_new)
summary(multi_intercept_only)

# With a p-value (Pr(Chi)) of 0.20 < 0.05, the change in deviance is not significant.
# Hence there is no improvement in the model with the variables.
anova(multi_intercept_only, mult_result, test='Chisq')
#########################################################################################################


### Total Hip Feature Selection #########################################################################
# run the RFE algorithm - Total Hip
set.seed(7)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
results_total_hip <- rfe(ds_new, ds$TH_status, sizes = seq_len(ncol(ds_new)), rfeControl=control)
# summarize the results - BMI, Age
print(results_total_hip)
# list the chosen features
predictors(results_total_hip)
# plot the results
plot(results_total_hip, type=c("g", "o"))

mult_result <- multinom(ds$TH_status ~ BMI + Age, maxit = 10000, data=ds_new)

summary(mult_result)
tab_model(mult_result, show.est=FALSE)

# Test whether the variables have an effect by testing the intercept only
multi_intercept_only <- multinom(ds$TH_status ~ 1, data=ds_new)
summary(multi_intercept_only)

# With a p-value (Pr(Chi)) of 0.84 > 0.05, the change in deviance is not significant.
# Hence there is no improvement in the model with the variables.
anova(multi_intercept_only, mult_result, test='Chisq')
#########################################################################################################