## IMPORTANT NOTE #######################################################################################
## Execute mb_data_filtering.R - From line 1 to 57 and 59 to 63 #########################################
#########################################################################################################

# Load necessary libraries
pacman::p_load(caret, nnet, sjPlot)

# Filter by Gender - get Males only
ds_new <- ds[ds$Gender == 0, ]
cat("Rows:", nrow(ds_new),", Columns:", ncol(ds_new),"\n")

Status_BMD <- ds_new$Status_BMD
Status_Fractures <- ds_new$Status_Fractures
LS_status <- ds_new$LS_status
FN_status <- ds_new$FN_status
TH_status <- ds_new$TH_status

# Setting explanatory variables only
additional_columns <- ds_new[, c("Age", "BMI", "Smoking", "Alcohol")]
ds_new <- cbind(additional_columns, ds_new[, 23:ncol(ds)])
cat("Rows:", nrow(ds_new),", Columns:", ncol(ds_new),"\n")

summary(ds_new)

# Drop Smoking column as all values are 0.
ds_new <- ds_new[-3]


### Find which attributes are highly correlated - multicolliniarity analysis
# calculate correlation matrix
correlationMatrix <- round(cor(ds_new, method='spearman'), 2)
# find attributes that are highly corrected
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
print(highlyCorrelated)
# Not hightly correlated @ 0.5 - 1-3, 36, 52, 54, 55, 61-63, 65-66, 74, 143, 163
# Get column name by column index colnames(ds_new[x])

######## Calculate Status BMD from manually selecting the variable with lower correlation
ds_new_status_bmd <- cbind(Status_BMD, ds_new)

mult_result <- multinom(ds_new_status_bmd$Status_BMD ~
                          Age + BMI + Alcohol + LDLsize + DHA + Gly + His + Tyr + Lactate +
                            Acetate + Acetone + Albumin + XXL.VLDL.TG + XL.HDL.FC + S.HDL.CE,
                        maxit = 1000, data = ds_new_status_bmd)
summary(mult_result)

## Reduced the variables to 7 according to the standard error vs the model
## coefficients (STD err has to be less than half of model co-off) and tab model
mult_result <- multinom(ds_new_status_bmd$Status_BMD ~ Age + LDLsize + His + Tyr + Acetate + Acetone + XL.HDL.FC,
                        maxit = 1000, data = ds_new_status_bmd)
summary(mult_result)
tab_model(mult_result, show.est=FALSE)

# Test whether the variables have an effect by testing the intercept only
multi_intercept_only <- multinom(ds_new_status_bmd$Status_BMD ~ 1, data=ds_new_status_bmd)
summary(multi_intercept_only)

# With a p-value (Pr(Chi)) of 0 < 0.05, the change in deviance is significant.
# Hence there is an improvement in the model with the variables.
anova(multi_intercept_only, mult_result, test='Chisq')

attach(ds_new_status_bmd)

# Prediction estimates
predicted <- predict(mult_result, ds_new_status_bmd$Status_BMD, type = "probs")
head(predicted)
predicted <- predict(mult_result, ds_new_status_bmd$Status_BMD, type = "class")
head(predicted)

detach(ds_new_status_bmd)

# Convert prediction to classifications
conttable <- with(ds_new_status_bmd,table(Status_BMD, predicted))
conttable

# Calculate the percentages of the predictions
percentCorrect <- round(diag(conttable)/rowSums(conttable),2)
finaltable <- cbind(conttable, percentCorrect)
finaltable
#########################################################################################################

######## Calculate Fracture Status from manually selecting the variable with lower correlation
ds_new_status_fractures <- cbind(Status_Fractures, ds_new)

lm_result <- glm(ds_new_status_fractures$Status_Fractures ~
                          Age + BMI + Alcohol + LDLsize + DHA + Gly + His + Tyr + Lactate +
                            Acetate + Acetone + Albumin + XXL.VLDL.TG + XL.HDL.FC + S.HDL.CE,
                        family = 'binomial', data = ds_new_status_fractures)
summary(lm_result)

## Reduced the variables to 2 for glm to converge
lm_result <- glm(ds_new_status_fractures$Status_Fractures ~ Gly, family = 'binomial', data = ds_new_status_fractures)
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

######## Calculate Lumbar Spine BMD from manually selecting the variable with lower correlation
ds_new_ls_bmd <- cbind(LS_status, ds_new)

mult_result <- multinom(LS_status ~
                          Age + BMI + Alcohol + LDLsize + DHA + Gly + His + Tyr + Lactate +
                            Acetate + Acetone + Albumin + XXL.VLDL.TG + XL.HDL.FC + S.HDL.CE,
                        maxit = 1000, data = ds_new_ls_bmd)
summary(mult_result)

## Reduced the variables to 6 according to the standard error vs the model
## coefficients (STD err has to be less than half of model co-off) and tab model
mult_result <- multinom(LS_status ~ LDLsize + Gly + His + Acetone + XL.HDL.FC + S.HDL.CE,
                        maxit = 1000, data = ds_new_ls_bmd)
summary(mult_result)
tab_model(mult_result, show.est=FALSE)

# Test whether the variables have an effect by testing the intercept only
multi_intercept_only <- multinom(ds_new_ls_bmd$LS_status ~ 1, data=ds_new_ls_bmd)
summary(multi_intercept_only)

# With a p-value (Pr(Chi)) of 0 < 0.05, the change in deviance is significant.
# Hence there is an improvement in the model with the variables.
anova(multi_intercept_only, mult_result, test='Chisq')

attach(ds_new_ls_bmd)

# Prediction estimates
predicted <- predict(mult_result, ds_new_ls_bmd$LS_status, type = "probs")
head(predicted)
predicted <- predict(mult_result, ds_new_ls_bmd$LS_status, type = "class")
head(predicted)

detach(ds_new_ls_bmd)

# Convert prediction to classifications
conttable <- with(ds_new_ls_bmd,table(LS_status, predicted))
conttable

# Calculate the percentages of the predictions
percentCorrect <- round(diag(conttable)/rowSums(conttable),2)
finaltable <- cbind(conttable, percentCorrect)
finaltable
#########################################################################################################

######## Calculate Femoral Neck BMD from manually selecting the variable with lower correlation
ds_new_fn_bmd <- cbind(FN_status, ds_new)

mult_result <- multinom(FN_status ~
                          Age + BMI + Alcohol + LDLsize + DHA + Gly + His + Tyr + Lactate +
                            Acetate + Acetone + Albumin + XXL.VLDL.TG + XL.HDL.FC + S.HDL.CE,
                        maxit = 1000, data = ds_new_fn_bmd)
summary(mult_result)

## Reduced the variables to 6 according to the standard error vs the model
## coefficients (STD err has to be less than half of model co-off) and tab model
mult_result <- multinom(FN_status ~  BMI + Alcohol + LDLsize + His + Acetone + XXL.VLDL.TG + XL.HDL.FC + S.HDL.CE,
                        maxit = 1000, data = ds_new_fn_bmd)
summary(mult_result)
tab_model(mult_result, show.est=FALSE)

# Test whether the variables have an effect by testing the intercept only
multi_intercept_only <- multinom(ds_new_fn_bmd$FN_status ~ 1, data=ds_new_fn_bmd)
summary(multi_intercept_only)

# With a p-value (Pr(Chi)) of 0 < 0.05, the change in deviance is significant.
# Hence there is an improvement in the model with the variables.
anova(multi_intercept_only, mult_result, test='Chisq')

attach(ds_new_fn_bmd)

# Prediction estimates
predicted <- predict(mult_result, ds_new_fn_bmd$FN_status, type = "probs")
head(predicted)
predicted <- predict(mult_result, ds_new_fn_bmd$FN_status, type = "class")
head(predicted)

detach(ds_new_fn_bmd)

# Convert prediction to classifications
conttable <- with(ds_new_fn_bmd,table(FN_status, predicted))
conttable

# Calculate the percentages of the predictions
percentCorrect <- round(diag(conttable)/rowSums(conttable),2)
finaltable <- cbind(conttable, percentCorrect)
finaltable
#########################################################################################################

######## Calculate Total Hip BMD from manually selecting the variable with lower correlation
ds_new_th_bmd <- cbind(TH_status, ds_new)

mult_result <- multinom(TH_status ~
                          Age + BMI + Alcohol + LDLsize + DHA + Gly + His + Tyr + Lactate +
                            Acetate + Acetone + Albumin + XXL.VLDL.TG + XL.HDL.FC + S.HDL.CE,
                        maxit = 1000, data = ds_new_th_bmd)
summary(mult_result)

## Reduced the variables to 6 according to the standard error vs the model
## coefficients (STD err has to be less than half of model co-off) and tab model
mult_result <- multinom(TH_status ~  Age + LDLsize + Gly + His + Tyr + Lactate + Acetone + XL.HDL.FC + S.HDL.CE,
                        maxit = 1000, data = ds_new_th_bmd)
summary(mult_result)
tab_model(mult_result, show.est=FALSE)

# Test whether the variables have an effect by testing the intercept only
multi_intercept_only <- multinom(ds_new_th_bmd$TH_status ~ 1, data=ds_new_th_bmd)
summary(multi_intercept_only)

# With a p-value (Pr(Chi)) of 0 < 0.05, the change in deviance is significant.
# Hence there is an improvement in the model with the variables.
anova(multi_intercept_only, mult_result, test='Chisq')

attach(ds_new_th_bmd)

# Prediction estimates
predicted <- predict(mult_result, ds_new_th_bmd$TH_status, type = "probs")
head(predicted)
predicted <- predict(mult_result, ds_new_th_bmd$TH_status, type = "class")
head(predicted)

detach(ds_new_th_bmd)

# Convert prediction to classifications
conttable <- with(ds_new_th_bmd,table(TH_status, predicted))
conttable

# Calculate the percentages of the predictions
percentCorrect <- round(diag(conttable)/rowSums(conttable),2)
finaltable <- cbind(conttable, percentCorrect)
finaltable
#########################################################################################################