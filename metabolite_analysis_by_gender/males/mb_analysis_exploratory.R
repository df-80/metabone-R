## IMPORTANT NOTE #######################################################################################
## Execute mb_data_filtering.R - From line 1 to 31 and line 42 to 65 ####################################
#########################################################################################################

# Load necessary libraries
pacman::p_load(ggplot2, reshape2)

## Define functions #####################################################################################
# define functions to get p-values for Kruskall-Wallis, Mann-Whitney and Pearson's
kwFunc <- function(var1, var2, data) {
  result <- kruskal.test(data[ ,var1], data[ ,var2], na.rm=TRUE)
  data.frame(var1, var2, result["p.value"], stringsAsFactors=FALSE)
}

mwFunc <- function(var1, var2, data) {
  result <- wilcox.test(data[ ,var1] ~ data[ ,var2], exact = FALSE)
  data.frame(var1, var2, result["p.value"], stringsAsFactors=FALSE)
}

corrFunc <- function(var1, var2, data, method) {
  result <- cor.test(data[ ,var1], data[ ,var2], na.rm=TRUE, method=method)
  data.frame(var1, var2, result[c("estimate","p.value","statistic","method")], stringsAsFactors=FALSE)
}

# Filter by Gender - get Males only
ds <- ds[ds$Gender == 'Male',]
cat("Rows:", nrow(ds),", Columns:", ncol(ds),"\n")

## Start Analysing data by applying MW or KW ############################################################
## Normality test of values between column 23 and column 193
## - p-value < 0.05 -> reject normality
## - p-value > 0.05 -> does not reject normality

# Apply Shapiro normality test on columns 23 to 193
lshap <- lapply(ds[23:ncol(ds)], shapiro.test)
# Extract relevant info from returned list of object
lres  <- sapply(lshap, `[`, c("statistic","p.value"))
# Transpose results list
shapiro_result <- t(lres)
write.csv(shapiro_result, "metabone-R/results/gender/males/shapiro.csv", row.names = TRUE)

## Pearson's Correlation Analysis
##  - p-value < 0.05 -> Correlation is significant
##  - p-value > 0.05 -> No correlation
#nd_ds <- ds[, c(1:22, 30:32, 34, 37, 64:65, 71:82, 84, 86:91, 94:102, 107:110, 113:114, 121, 149, 156:157, 163, 170, 177, 179:180, 183)]
nd_ds <- ds[, c(1:22, 23:29, 33, 35:36, 38:63, 66:70, 83, 85, 92:93, 103:106, 111:112, 115:120, 122:148, 150:155, 158:162, 164:169, 171:176, 178, 181:182, 184:191)]
vars <- data.frame(v1=names(nd_ds[8]), v2=names(nd_ds[,23:ncol(nd_ds)]))
tmp_ds <- nd_ds[,c(8, 23:ncol(nd_ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(corrFunc, vars[,1], vars[,2], MoreArgs=list(data=tmp_ds, method="pearson"), SIMPLIFY=FALSE))
pearson_result <- corrs[2:5]
write.csv(pearson_result, "metabone-R/results/gender/males/pearsonFractureStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# corrs[4] <- p.adjust(corrs$p.value, method = "bonferroni", n = length(corrs$p.value))
# pearson_result <- corrs[2:5]
# write.csv(pearson_result, "metabone-R/results/gender/males/pearsonFractureStatusvsMetabolitesAdjusted.csv", row.names = TRUE)

## Pearson's Correlation Analysis on all vars
##  - p-value < 0.05 -> Correlation is significant
##  - p-value > 0.05 -> No correlation
nd_ds <- ds
vars <- data.frame(v1=names(nd_ds[8]), v2=names(nd_ds[,23:ncol(nd_ds)]))
tmp_ds <- nd_ds[,c(8, 23:ncol(nd_ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(corrFunc, vars[,1], vars[,2], MoreArgs=list(data=tmp_ds, method="pearson"), SIMPLIFY=FALSE))
pearson_result <- corrs[2:5]
write.csv(pearson_result, "metabone-R/results/gender/males/pearsonFractureStatusvsMetabolites2.csv", row.names = TRUE)

## Mann-Whitney/Kruskall Wallis Analysis
## - p-value < 0.05 -> Correlation is significant
## - p-value > 0.05 -> No correlation

## Apply Mann-Whitney on Fracture Status [8]
vars <- data.frame(v1=names(ds[8]), v2=names(ds[,23:ncol(ds)]))
tmp_ds <- ds[,c(8, 23:ncol(ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(mwFunc, vars[,2], vars[,1], MoreArgs=list(data=tmp_ds), SIMPLIFY=FALSE))
mann_whitney_result <- corrs[c(1,3)]
write.csv(mann_whitney_result, "metabone-R/results/gender/males/mw/FractureStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# mann_whitney_result$p.value <- p.adjust(mann_whitney_result$p.value, method = "bonferroni", n = length(mann_whitney_result$p.value))
# write.csv(mann_whitney_result, "metabone-R/results/gender/males/mw/FractureStatusvsMetabolitesAdjusted.csv", row.names = TRUE)

## Apply Kruskall-Wallis on Status_BMD [7]
vars <- data.frame(v1=names(ds[7]), v2=names(ds[,23:ncol(ds)]))
tmp_ds <- ds[,c(7, 23:ncol(ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(kwFunc, vars[,2], vars[,1], MoreArgs=list(data=tmp_ds), SIMPLIFY=FALSE))
kruskal_wallis_result <- corrs[2:3]
write.csv(kruskal_wallis_result, "metabone-R/results/gender/males/kw/BMDStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# kruskal_wallis_result$p.value <- p.adjust(kruskal_wallis_result$p.value, method = "bonferroni", n = length(kruskal_wallis_result$p.value))
# write.csv(kruskal_wallis_result, "metabone-R/results/gender/males/kw/BMDStatusvsMetabolitesAdjusted.csv", row.names = TRUE)

## Apply Kruskall-Wallis on LS Status [10]
vars <- data.frame(v1=names(ds[10]), v2=names(ds[,23:ncol(ds)]))
tmp_ds <- ds[,c(10, 23:ncol(ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(kwFunc, vars[,2], vars[,1], MoreArgs=list(data=tmp_ds), SIMPLIFY=FALSE))
kruskal_wallis_result <- corrs[c(2:3)]
write.csv(kruskal_wallis_result, "metabone-R/results/gender/males/kw/LsStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# kruskal_wallis_result$p.value <- p.adjust(kruskal_wallis_result$p.value, method = "bonferroni", n = length(kruskal_wallis_result$p.value))
# write.csv(kruskal_wallis_result, "metabone-R/results/gender/males/kw/LsStatusvsMetabolitesAdjusted.csv", row.names = TRUE)

## Apply Kruskall-Wallis on FN Status [12]
vars <- data.frame(v1=names(ds[12]), v2=names(ds[,23:ncol(ds)]))
tmp_ds <- ds[,c(12, 23:ncol(ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(kwFunc, vars[,2], vars[,1], MoreArgs=list(data=tmp_ds), SIMPLIFY=FALSE))
kruskal_wallis_result <- corrs[2:3]
write.csv(kruskal_wallis_result, "metabone-R/results/gender/males/kw/FnStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# kruskal_wallis_result$p.value <- p.adjust(kruskal_wallis_result$p.value, method = "bonferroni", n = length(kruskal_wallis_result$p.value))
# write.csv(kruskal_wallis_result, "metabone-R/results/gender/males/kw/FnStatusvsMetabolitesAdjusted.csv", row.names = TRUE)

## Apply Kruskall-Wallis on TH Status [14]
vars <- data.frame(v1=names(ds[14]), v2=names(ds[,23:ncol(ds)]))
tmp_ds <- ds[,c(14, 23:ncol(ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(kwFunc, vars[,2], vars[,1], MoreArgs=list(data=tmp_ds), SIMPLIFY=FALSE))
kruskal_wallis_result <- corrs[2:3]
write.csv(kruskal_wallis_result, "metabone-R/results/gender/males/kw/ThStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# kruskal_wallis_result$p.value <- p.adjust(kruskal_wallis_result$p.value, method = "bonferroni", n = length(kruskal_wallis_result$p.value))
# write.csv(kruskal_wallis_result, "metabone-R/results/gender/males/kw/ThStatusvsMetabolitesAdjusted.csv", row.names = TRUE)


## Plot Boxplots ##########################################################################################
## Status Fractures #######################################################################################

tmp_ds <- ds[,c(8, 23:82)]
ds.m <- melt(tmp_ds, id.var = "Status_Fractures")
ggplot(data = ds.m, aes(x=Status_Fractures, y=value)) +
  geom_boxplot(aes(fill = Status_Fractures)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Status Fractures vs Metabolites - 1")

tmp_ds <- ds[,c(8, 83:142)]
ds.m <- melt(tmp_ds, id.var = "Status_Fractures")
ggplot(data = ds.m, aes(x=Status_Fractures, y=value)) +
  geom_boxplot(aes(fill = Status_Fractures)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Status Fractures vs Metabolites - 2")

tmp_ds <- ds[,c(8, 143:ncol(ds))]
ds.m <- melt(tmp_ds, id.var = "Status_Fractures")
ggplot(data = ds.m, aes(x=Status_Fractures, y=value)) +
  geom_boxplot(aes(fill = Status_Fractures)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Status Fractures vs Metabolites - 3")

## Most significant metabolites vs Status Fractures (without outliers) ####################################
## Status Fractures vs Gly ################################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = Gly, na.rm = TRUE, fill = Status_Fractures)) +
#  geom_boxplot( , na.rm = TRUE) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Gly, c(0, 0.9) * 1.1, na.rm = TRUE)) +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs Gly", x = "Status_Fractures", y = "Gly")
p

kruskal.test(ds$Gly, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(Gly ~ Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs Total.CE ###########################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = Total.CE, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Total.CE, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs Total.CE", x = "Status_Fractures", y = "Total.CE")
p

kruskal.test(ds$Total.CE, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(Total.CE ~ ds$Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs Sphingomyelins #####################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = Sphingomyelins, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Sphingomyelins, c(0, 0.9) * 1.05, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "Sphingomyelins") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs Sphingomyelins", x = "Status_Fractures", y = "Sphingomyelins")
p

kruskal.test(ds$Sphingomyelins, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(Sphingomyelins ~ Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs Glucose ############################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = Glucose, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Glucose, c(0, 0.9) * 1.09, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "Glucose") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs Glucose", x = "Status_Fractures", y = "Glucose")
p

kruskal.test(ds$Glucose, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(Glucose ~ Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs Total.C ############################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = Total.C, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Total.C, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "Total.C") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs Total.C", x = "Status_Fractures", y = "Total.C")
p

kruskal.test(ds$Total.C, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(Total.C ~ Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs Total.FC ############################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = Total.FC, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Total.FC, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "Total.FC") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs Total.FC", x = "Status_Fractures", y = "Total.FC")
p

kruskal.test(ds$Total.FC, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(Total.FC ~ Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs ApoB ###############################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = ApoB, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$ApoB, c(0, 0.9) * 1.1, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "ApoB") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs ApoB", x = "Status_Fractures", y = "ApoB")
p

kruskal.test(ds$ApoB, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(ApoB ~ Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs Total.L ############################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = Total.L, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Total.L, c(0, 0.9) * 1.1, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "Total.L") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs Total.L", x = "Status_Fractures", y = "Total.L")
p

kruskal.test(ds$Total.L, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(Total.L ~ Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs Omega.6 ############################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = Omega.6, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Omega.6, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "Omega.6") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs Omega.6", x = "Status_Fractures", y = "Omega.6")
p

kruskal.test(ds$Omega.6, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(Omega.6 ~ Status_Fractures, data = ds, exact = FALSE)


## Status BMD #############################################################################################
tmp_ds <- ds[,c(7, 23:82)]
ds.m <- melt(tmp_ds, id.var = "Status_BMD")
ggplot(data = ds.m, aes(x=Status_BMD, y=value)) +
  geom_boxplot(aes(fill = Status_BMD)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Status BMD vs Metabolites - 1")

tmp_ds <- ds[,c(7, 83:142)]
ds.m <- melt(tmp_ds, id.var = "Status_BMD")
ggplot(data = ds.m, aes(x=Status_BMD, y=value)) +
  geom_boxplot(aes(fill = Status_BMD)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Status BMD vs Metabolites - 2")

tmp_ds <- ds[,c(7, 143:ncol(ds))]
ds.m <- melt(tmp_ds, id.var = "Status_BMD")
ggplot(data = ds.m, aes(x=Status_BMD, y=value)) +
  geom_boxplot(aes(fill = Status_BMD)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Status BMD vs Metabolites - 3")

## Most significant metabolites vs Status BMD (without outliers) ##########################################
## Status BMD vs Sphingomyelins ###########################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = Sphingomyelins, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Sphingomyelins, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "Sphingomyelins") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs Sphingomyelins", x = "Status_BMD", y = "Sphingomyelins")
p

by(ds$Sphingomyelins, ds$Status_BMD, shapiro.test)
kruskal.test(ds$Sphingomyelins, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$Sphingomyelins, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## Status BMD vs Gly ######################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = Gly, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Gly, c(0, 0.9) * 1.1, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "Gly") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs Gly", x = "Status_BMD", y = "Gly")
p

by(ds$Gly, ds$Status_BMD, shapiro.test)
kruskal.test(ds$Gly, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$Gly, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## Status BMD vs Total.C ##################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = Total.C, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Total.C, c(0, 0.9) * 1.05, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "Total.C") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs Total.C", x = "Status_BMD", y = "Total.C")
p

by(ds$Total.C, ds$Status_BMD, shapiro.test)
kruskal.test(ds$Total.C, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$Total.C, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## Status BMD vs ApoB #####################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = ApoB, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$ApoB, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "ApoB") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs ApoB", x = "Status_BMD", y = "ApoB")
p

by(ds$ApoB, ds$Status_BMD, shapiro.test)
kruskal.test(ds$ApoB, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$ApoB, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## Status BMD vs Glycerol ################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = Glycerol, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Glycerol, c(0, 0.9) * 1.10, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "Glycerol") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs Glycerol", x = "Status_BMD", y = "Glycerol")
p

by(ds$Glycerol, ds$Status_BMD, shapiro.test)
kruskal.test(ds$Glycerol, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$Glycerol, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## Status BMD vs Glucose ##################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = Glucose, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Glucose, c(0, 0.9) * 1.05, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "Glucose") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs Glucose", x = "Status_BMD", y = "Glucose")
p

by(ds$Glucose, ds$Status_BMD, shapiro.test)
kruskal.test(ds$Glucose, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$Glucose, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## Status BMD vs His ######################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = His, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$His, c(0, 0.9) * 1.06, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "His") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs His", x = "Status_BMD", y = "His")
p

by(ds$His, ds$Status_BMD, shapiro.test)
kruskal.test(ds$His, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$His, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## Status BMD vs Tyr ######################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = Leu, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Leu, c(0, 0.9) * 1.10, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "Leu") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs Leu", x = "Status_BMD", y = "Leu")
p

by(ds$Leu, ds$Status_BMD, shapiro.test)
kruskal.test(ds$Leu, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$Leu, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")


## Lumbar Spine ###########################################################################################

tmp_ds <- ds[,c(10, 23:82)]
ds.m <- melt(tmp_ds, id.var = "LS_status")
ggplot(data = ds.m, aes(x=LS_status, y=value)) +
  geom_boxplot(aes(fill = LS_status)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Lumbar Spine vs Metabolites - 1")

tmp_ds <- ds[,c(10, 83:142)]
ds.m <- melt(tmp_ds, id.var = "LS_status")
ggplot(data = ds.m, aes(x=LS_status, y=value)) +
  geom_boxplot(aes(fill = LS_status)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Lumbar Spine vs Metabolites - 2")

tmp_ds <- ds[,c(10, 143:ncol(ds))]
ds.m <- melt(tmp_ds, id.var = "LS_status")
ggplot(data = ds.m, aes(x=LS_status, y=value)) +
  geom_boxplot(aes(fill = LS_status)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Lumbar Spine vs Metabolites - 3")

## Most significant metabolites vs LS Status (without outliers) ##########################################
## LS status vs Glucose ##################################################################################
p <- ggplot(ds, aes(x = LS_status, y = Glucose, na.rm = TRUE, fill = LS_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Glucose, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "LS_status", y = "Glucose") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "LS_status vs Glucose", x = "LS_status", y = "Glucose")
p

by(ds$Glucose, ds$LS_status, shapiro.test)
kruskal.test(ds$Glucose, ds$LS_status, rm.na=TRUE)
pairwise.wilcox.test(ds$Glucose, ds$LS_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

p <- ggplot(ds, aes(x = LS_status, y = His, na.rm = TRUE, fill = LS_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$His, c(0, 0.9) * 1.08, na.rm = TRUE)) +
  labs(x = "LS_status", y = "His") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "LS_status vs His", x = "LS_status", y = "His")
p

by(ds$His, ds$LS_status, shapiro.test)
kruskal.test(ds$His, ds$LS_status, rm.na=TRUE)
pairwise.wilcox.test(ds$His, ds$LS_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

p <- ggplot(ds, aes(x = LS_status, y = Ile, na.rm = TRUE, fill = LS_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Ile, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "LS_status", y = "Ile") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "LS_status vs Ile", x = "LS_status", y = "Ile")
p

by(ds$Ile, ds$LS_status, shapiro.test)
kruskal.test(ds$Ile, ds$LS_status, rm.na=TRUE)
pairwise.wilcox.test(ds$Ile, ds$LS_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")


## Total Hip ##############################################################################################

tmp_ds <- ds[,c(14, 23:82)]
ds.m <- melt(tmp_ds, id.var = "TH_status")
ggplot(data = ds.m, aes(x=TH_status, y=value)) +
  geom_boxplot(aes(fill = TH_status)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Total Hip vs Metabolites - 1")

tmp_ds <- ds[,c(14, 83:142)]
ds.m <- melt(tmp_ds, id.var = "TH_status")
ggplot(data = ds.m, aes(x=TH_status, y=value)) +
  geom_boxplot(aes(fill = TH_status)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Total Hip vs Metabolites - 2")

tmp_ds <- ds[,c(14, 143:ncol(ds))]
ds.m <- melt(tmp_ds, id.var = "TH_status")
ggplot(data = ds.m, aes(x=TH_status, y=value)) +
  geom_boxplot(aes(fill = TH_status)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Total Hip vs Metabolites - 3")

## Most significant metabolites vs TH Status (without outliers) ###########################################
## TH status vs Gly #######################################################################################
p <- ggplot(ds, aes(x = TH_status, y = Gly, na.rm = TRUE, fill = TH_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Gly, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "TH_status", y = "Gly") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "TH_status vs Gly", x = "TH_status", y = "Gly")
p

by(ds$Gly, ds$TH_status, shapiro.test)
kruskal.test(ds$Gly, ds$TH_status, rm.na=TRUE)
pairwise.wilcox.test(ds$Gly, ds$TH_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## TH status vs Acetoacetate ##############################################################################
p <- ggplot(ds, aes(x = TH_status, y = Acetoacetate, na.rm = TRUE, fill = TH_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Acetoacetate, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "TH_status", y = "Acetoacetate") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "TH_status vs Acetoacetate", x = "TH_status", y = "Acetoacetate")
p

by(ds$Acetoacetate, ds$TH_status, shapiro.test)
kruskal.test(ds$Acetoacetate, ds$TH_status, rm.na=TRUE)
pairwise.wilcox.test(ds$Acetoacetate, ds$TH_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")


## Femoral Neck ###########################################################################################

tmp_ds <- ds[,c(12, 23:82)]
ds.m <- melt(tmp_ds, id.var = "FN_status")
ggplot(data = ds.m, aes(x=FN_status, y=value)) +
  geom_boxplot(aes(fill = FN_status)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Femoral Neck vs Metabolites - 1")

tmp_ds <- ds[,c(12, 83:142)]
ds.m <- melt(tmp_ds, id.var = "FN_status")
ggplot(data = ds.m, aes(x=FN_status, y=value)) +
  geom_boxplot(aes(fill = FN_status)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Femoral Neck vs Metabolites - 2")

tmp_ds <- ds[,c(12, 143:ncol(ds))]
ds.m <- melt(tmp_ds, id.var = "FN_status")
ggplot(data = ds.m, aes(x=FN_status, y=value)) +
  geom_boxplot(aes(fill = FN_status)) + facet_wrap( ~ variable, ncol = 8, scales="free") + ggtitle("Femoral Neck vs Metabolites - 3")

## Most significant metabolites vs FN Status (without outliers) ########################################
## FN status vs Sphingomyelins #########################################################################
p <- ggplot(ds, aes(x = FN_status, y = Sphingomyelins, na.rm = TRUE, fill = FN_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Sphingomyelins, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "FN_status", y = "Sphingomyelins") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "FN_status vs Sphingomyelins", x = "FN_status", y = "Sphingomyelins")
p

by(ds$Sphingomyelins, ds$FN_status, shapiro.test)
kruskal.test(ds$Sphingomyelins, ds$FN_status, rm.na=TRUE)
pairwise.wilcox.test(ds$Sphingomyelins, ds$FN_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## FN status vs SFA #######################################################################################
p <- ggplot(ds, aes(x = FN_status, y = SFA, na.rm = TRUE, fill = FN_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$SFA, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "FN_status", y = "SFA") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "FN_status vs SFA", x = "FN_status", y = "SFA")
p

by(ds$SFA, ds$FN_status, shapiro.test)
kruskal.test(ds$SFA, ds$FN_status, rm.na=TRUE)
pairwise.wilcox.test(ds$SFA, ds$FN_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")