# Males data is not sufficient for SVM.
# ## IMPORTANT NOTE #######################################################################################
# ## Execute mb_data_filtering.R - From line 1 to 32 and 44 to 51 #########################################
# #########################################################################################################
#
# # Load necessary libraries
# pacman::p_load(caret, e1071, nnet, ggplot2, sjPlot)
#
# ds <- data.frame(dataset['Age'], dataset['BMI'], ds)
#
# ## Scale data - Use caret library to MinMax scale all columns
# process <- preProcess(as.data.frame(ds), method = "range")
# norm_scale <- predict(process, as.data.frame(ds))
#
# ###################################################################################################
# ## Status Fractures ###############################################################################
# ###################################################################################################
#
# # Add Status_Fractures, Gender, Smoking, Alcohol, Age and BMI to a new dataframe
# ds_new <- data.frame(
#   dataset['Status_Fractures'], dataset['Gender'], dataset['Smoking'],  dataset['Alcohol'], norm_scale['Age'], norm_scale['BMI'],
#   norm_scale['Gly'], norm_scale['L.LDL.FC'], norm_scale['LDL.FC'], norm_scale['M.LDL.FC'], norm_scale['Clinical.LDL.C'],
#   norm_scale['Total.CE'], norm_scale['L.LDL.PL'], norm_scale['L.LDL.C'], norm_scale['L.LDL.CE'],
#   norm_scale['S.LDL.FC'], norm_scale['Sphingomyelins'], norm_scale['LDL.PL'], norm_scale['L.LDL.L'], norm_scale['LDL.C'],
#   norm_scale['Glucose'], norm_scale['S.LDL.C'], norm_scale['M.VLDL.CE'], norm_scale['S.LDL.P'], norm_scale['Total.C'],
#   norm_scale['LDL.CE'], norm_scale['LDL.L'], norm_scale['IDL.C'], norm_scale['IDL.FC'], norm_scale['L.LDL.P'],
#   norm_scale['M.LDL.C'], norm_scale['S.LDL.L'], norm_scale['S.LDL.CE'], norm_scale['Total.FC'], norm_scale['S.LDL.PL'],
#   norm_scale['IDL.CE'], norm_scale['M.LDL.PL'], norm_scale['non.HDL.C'], norm_scale['M.VLDL.C'], norm_scale['M.LDL.L'],
#   norm_scale['M.LDL.CE'], norm_scale['IDL.L'], norm_scale['IDL.PL'], norm_scale['LDL.P'], norm_scale['IDL.P'],
#   norm_scale['ApoB'], norm_scale['Remnant.C'], norm_scale['M.HDL.TG'], norm_scale['M.VLDL.FC'], norm_scale['S.VLDL.FC'],
#   norm_scale['Total.L'], norm_scale['Omega.6'], norm_scale['HDL.TG'], norm_scale['M.LDL.P'])
#
# # Filter by Gender - get Males only
# ds_new <- ds_new[ds_new$Gender == 0, ]
#
# # Drop Smoking and Gender columns as all values are 0.
# ds_new <- ds_new[c(-2, -3)]
#
# # Change Status Fractures to categorical variable
# ds_new$Status_Fractures <- factor(ds_new$Status_Fractures, levels=c(0,1), labels=c("No", "Yes"))
#
# # Remove columns which have more than 25 NAs
# ds_new <- ds_new[, colSums(is.na(ds_new)) <= 24]
# ds.noNA <- ds_new %>% na.omit()
#
# summary(ds.noNA)
# balanced.ds <- ds.noNA %>% group_by(Status_Fractures) %>% slice_sample(n = 5)
#
# ## 70% of the sample size
# train_size <- floor(0.7 * nrow(balanced.ds))
#
# set.seed(55)
# in_rows <- sample(seq_len(nrow(balanced.ds)), size = train_size, replace = FALSE)
#
# train <- balanced.ds[in_rows, ]
# test  <- balanced.ds[-in_rows, ]
# attach(train)
#
# mult_result <- glm(Status_Fractures ~
#                      Age + BMI + Gly + L.LDL.FC + LDL.FC + M.LDL.FC + Clinical.LDL.C + Total.CE + L.LDL.PL + L.LDL.C + L.LDL.CE +
#                        S.LDL.FC + Sphingomyelins + LDL.PL + L.LDL.L + LDL.C + Glucose + S.LDL.C + M.VLDL.CE + S.LDL.P + Total.C +
#                        LDL.CE + LDL.L + IDL.C + IDL.FC + L.LDL.P + M.LDL.C + S.LDL.L + S.LDL.CE + Total.FC + S.LDL.PL + IDL.CE +
#                        M.LDL.PL + non.HDL.C + M.VLDL.C + M.LDL.L + M.LDL.CE + IDL.L + IDL.PL + LDL.P + IDL.P + ApoB + Remnant.C +
#                        M.HDL.TG + M.VLDL.FC + S.VLDL.FC + Total.L + Omega.6 + HDL.TG + M.LDL.P,
#                    family = 'binomial', data=train)
# mult_result
#
# predict <- predict(mult_result, newdata=test, type="response")
# predict.cat <- cut(predict,breaks=c(0, 0.5, 1), labels=c("No fractures", "Fractures"))
# summary(predict.cat)
# summary(test$Status_Fractures)
# table(predict.cat, test$Status_Fractures)
#
# svmfit  <- svm(Status_Fractures ~
#                  Age + BMI + Gly + L.LDL.FC + LDL.FC + M.LDL.FC + Clinical.LDL.C + Total.CE + L.LDL.PL + L.LDL.C + L.LDL.CE +
#                    S.LDL.FC + Sphingomyelins + LDL.PL + L.LDL.L + LDL.C + Glucose + S.LDL.C + M.VLDL.CE + S.LDL.P + Total.C +
#                    LDL.CE + LDL.L + IDL.C + IDL.FC + L.LDL.P + M.LDL.C + S.LDL.L + S.LDL.CE + Total.FC + S.LDL.PL + IDL.CE +
#                    M.LDL.PL + non.HDL.C + M.VLDL.C + M.LDL.L + M.LDL.CE + IDL.L + IDL.PL + LDL.P + IDL.P + ApoB + Remnant.C +
#                    M.HDL.TG + M.VLDL.FC + S.VLDL.FC + Total.L + Omega.6 + HDL.TG + M.LDL.P,
#                data = train, kernel = "linear", cost = 10, scale = FALSE)
#
# predict <- predict(svmfit, newdata=test, type="response")
# actual  <- test$Status_Fractures
# confusionMatrix(table(predict, actual))
#
#
# ###################################################################################################
# ## Status BMD #####################################################################################
# ###################################################################################################
#
# # Add Status_Fractures, Gender, Smoking, Alcohol, Age and BMI to a new dataframe
# ds_new <- data.frame(
#   dataset['Status_BMD'], dataset['Gender'], dataset['Smoking'],  dataset['Alcohol'], norm_scale['Age'], norm_scale['BMI'],
#   norm_scale['Sphingomyelins'], norm_scale['Gly'], norm_scale['M.VLDL.FC'], norm_scale['S.LDL.P'], norm_scale['Total.C'],
#   norm_scale['S.LDL.CE'], norm_scale['LDL.P'], norm_scale['M.VLDL.CE'], norm_scale['M.VLDL.C'], norm_scale['L.LDL.P'],
#   norm_scale['ApoB'], norm_scale['S.LDL.L'], norm_scale['S.LDL.PL'], norm_scale['Glycerol'], norm_scale['Glucose'],
#   norm_scale['S.LDL.C'], norm_scale['non.HDL.C'], norm_scale['L.LDL.CE'], norm_scale['IDL.P'], norm_scale['S.VLDL.FC'],
#   norm_scale['L.LDL.C'], norm_scale['M.VLDL.PL'], norm_scale['L.LDL.FC'], norm_scale['Total.CE'], norm_scale['Clinical.LDL.C'],
#   norm_scale['HDLsize'], norm_scale['M.LDL.P'], norm_scale['His'], norm_scale['Total.L'], norm_scale['L.LDL.L'],
#   norm_scale['IDL.FC'], norm_scale['LDL.CE'], norm_scale['LDL.L'], norm_scale['LDL.C'], norm_scale['M.LDL.CE'],
#   norm_scale['LDL.PL'], norm_scale['L.LDL.PL'], norm_scale['M.LDL.L'], norm_scale['M.VLDL.P'])
#
# # Filter by Gender - get Males only
# ds_new <- ds_new[ds_new$Gender == 0, ]
#
# # Drop Smoking and Gender columns as all values are 0.
# ds_new <- ds_new[c(-2, -3)]
#
# # Change BMD to categorical variable
# ds_new$Status_BMD <- factor(ds_new$Status_BMD, levels=c(0,1,2), labels=c("Normal", "Osteopenic", "Osteoporotic"))
#
# # Remove columns which have more than 25 NAs
# ds.noNA <- ds_new %>% na.omit()
#
# summary(ds.noNA)
# balanced.ds <- ds.noNA %>% group_by(Status_BMD) %>% slice_sample(n = 5)
#
# ## 70% of the sample size
# train_size <- floor(0.7 * nrow(balanced.ds))
#
# set.seed(33)
# in_rows <- sample(seq_len(nrow(balanced.ds)), size = train_size, replace = FALSE)
#
# train <- balanced.ds[in_rows, ]
# test  <- balanced.ds[-in_rows, ]
# attach(train)
#
# ## Multinomial Regression #############################################################################
# mult_result <- multinom(
#   train$Status_BMD ~ Age + BMI + Sphingomyelins + Gly + M.VLDL.FC + S.LDL.P + Total.C + S.LDL.CE + LDL.P + M.VLDL.CE + M.VLDL.C +
#     L.LDL.P + ApoB + S.LDL.L + S.LDL.PL + Glycerol + Glucose + S.LDL.C + non.HDL.C + L.LDL.CE + IDL.P + S.VLDL.FC + L.LDL.C +
#     M.VLDL.PL + L.LDL.FC + Total.CE + Clinical.LDL.C + HDLsize + M.LDL.P + His + Total.L + L.LDL.L + IDL.FC + LDL.CE + LDL.L +
#     LDL.C + M.LDL.CE + LDL.PL + L.LDL.PL + M.LDL.L + M.VLDL.P,
#   maxit = 100000, data=train)
#
# summary(mult_result)
# tab_model(mult_result, show.est=FALSE)
#
# predict <- predict(mult_result, newdata=test, type="class")
# actual <- test$Status_BMD
# confusionMatrix(table(predict, actual))
#
# ## SVM ################################################################################################
# svmfit  <- svm(
#   train$Status_BMD ~ Age + BMI + Sphingomyelins + Gly + M.VLDL.FC + S.LDL.P + Total.C + S.LDL.CE + LDL.P + M.VLDL.CE + M.VLDL.C +
#     L.LDL.P + ApoB + S.LDL.L + S.LDL.PL + Glycerol + Glucose + S.LDL.C + non.HDL.C + L.LDL.CE + IDL.P + S.VLDL.FC + L.LDL.C +
#     M.VLDL.PL + L.LDL.FC + Total.CE + Clinical.LDL.C + HDLsize + M.LDL.P + His + Total.L + L.LDL.L + IDL.FC + LDL.CE + LDL.L +
#     LDL.C + M.LDL.CE + LDL.PL + L.LDL.PL + M.LDL.L + M.VLDL.P,
#   data = train, kernel='radial', cost = 10, scale = TRUE)
#
# predict <- predict(svmfit, newdata=test, type="response")
# actual  <- test$Status_BMD
# confusionMatrix(table(predict, actual))
#
# ###################################################################################################
# ## Lumbar Spine ###################################################################################
# ###################################################################################################
#
# # # Add Status_Fractures, Gender, Smoking, Alcohol, Age and BMI to a new dataframe
# # ds_new <- data.frame(
# #   dataset['LS_status'], dataset['Gender'], dataset['Smoking'],  dataset['Alcohol'], norm_scale['Age'], norm_scale['BMI'],
# #   norm_scale['LDLsize'], norm_scale['Citrate'], norm_scale['His'])
# #
# # # Change Lumbar Spine status to categorical variable
# # ds_new$LS_status <- factor(ds_new$LS_status, levels=c(0,1,2), labels=c("Normal", "Osteopenic", "Osteoporotic"))
# #
# # # Remove columns which have more than 25 NAs
# # ds.noNA <- ds_new %>% na.omit()
# #
# # summary(ds.noNA)
# # balanced.ds <- ds.noNA %>% group_by(LS_status) %>% slice_sample(n = 44)
# #
# # ## 70% of the sample size
# # train_size <- floor(0.7 * nrow(balanced.ds))
# #
# # set.seed(47)
# # in_rows <- sample(seq_len(nrow(balanced.ds)), size = train_size, replace = FALSE)
# #
# # train <- balanced.ds[in_rows, ]
# # test  <- balanced.ds[-in_rows, ]
# # attach(train)
# #
# # ## Multinomial Regression #############################################################################
# # mult_result <- multinom(
# #   train$LS_status ~ Citrate + His,
# #   maxit = 100000, data=train)
# #
# # summary(mult_result)
# # tab_model(mult_result, show.est=FALSE)
# #
# # predict <- predict(mult_result, newdata=test, type="class")
# # actual <- test$LS_status
# # confusionMatrix(table(predict, actual))
# #
# # ## SVM ################################################################################################
# # svmfit  <- svm(
# #   train$LS_status ~
# #     Citrate + His, data = train, kernel='radial', cost = 10, scale = TRUE)
# # predict <- predict(svmfit, newdata=test, type="response")
# # actual  <- test$LS_status
# # confusionMatrix(table(predict, actual))