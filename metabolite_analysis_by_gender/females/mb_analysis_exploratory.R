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

# Filter by Gender - get Females only
ds <- ds[ds$Gender == 'Female',]
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
write.csv(shapiro_result, "metabone-R/results/gender/females/shapiro.csv", row.names = TRUE)

## Pearson's Correlation Analysis
##  - p-value < 0.05 -> Correlation is significant
##  - p-value > 0.05 -> No correlation
nd_ds <- ds[, c(1:22, 25, 34, 36, 41, 42, 44, 48, 50, 53, 54, 57, 58, 61, 62, 63, 93, 136:141, 143, 145:147, 152, 158, 160:162, 185:190)]
vars <- data.frame(v1=names(nd_ds[8]), v2=names(nd_ds[,23:ncol(nd_ds)]))
tmp_ds <- nd_ds[,c(8, 23:ncol(nd_ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(corrFunc, vars[,1], vars[,2], MoreArgs=list(data=tmp_ds, method="pearson"), SIMPLIFY=FALSE))
pearson_result <- corrs[2:5]
write.csv(pearson_result, "metabone-R/results/gender/females/pearsonFractureStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# corrs[4] <- p.adjust(corrs$p.value, method = "bonferroni", n = length(corrs$p.value))
# pearson_result <- corrs[2:5]
# write.csv(pearson_result, "metabone-R/results/mw_kw/pearsonFractureStatusvsMetabolitesAdjusted.csv", row.names = TRUE)

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
write.csv(pearson_result, "metabone-R/results/gender/females/pearsonFractureStatusvsMetabolites2.csv", row.names = TRUE)

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
write.csv(mann_whitney_result, "metabone-R/results/gender/females/mw/FractureStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# mann_whitney_result$p.value <- p.adjust(mann_whitney_result$p.value, method = "bonferroni", n = length(mann_whitney_result$p.value))
# write.csv(mann_whitney_result, "metabone-R/results/gender/females/mw/FractureStatusvsMetabolitesAdjusted.csv", row.names = TRUE)

## Apply Kruskall-Wallis on Status_BMD [7]
vars <- data.frame(v1=names(ds[7]), v2=names(ds[,23:ncol(ds)]))
tmp_ds <- ds[,c(7, 23:ncol(ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(kwFunc, vars[,2], vars[,1], MoreArgs=list(data=tmp_ds), SIMPLIFY=FALSE))
kruskal_wallis_result <- corrs[2:3]
write.csv(kruskal_wallis_result, "metabone-R/results/gender/females/kw/BMDStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# kruskal_wallis_result$p.value <- p.adjust(kruskal_wallis_result$p.value, method = "bonferroni", n = length(kruskal_wallis_result$p.value))
# write.csv(kruskal_wallis_result, "metabone-R/results/gender/females/kw/BMDStatusvsMetabolitesAdjusted.csv", row.names = TRUE)

## Apply Kruskall-Wallis on LS Status [10]
vars <- data.frame(v1=names(ds[10]), v2=names(ds[,23:ncol(ds)]))
tmp_ds <- ds[,c(10, 23:ncol(ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(kwFunc, vars[,2], vars[,1], MoreArgs=list(data=tmp_ds), SIMPLIFY=FALSE))
kruskal_wallis_result <- corrs[c(2:3)]
write.csv(kruskal_wallis_result, "metabone-R/results/gender/females/kw/LsStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# kruskal_wallis_result$p.value <- p.adjust(kruskal_wallis_result$p.value, method = "bonferroni", n = length(kruskal_wallis_result$p.value))
# write.csv(kruskal_wallis_result, "metabone-R/results/gender/females/kw/LsStatusvsMetabolitesAdjusted.csv", row.names = TRUE)

## Apply Kruskall-Wallis on FN Status [12]
vars <- data.frame(v1=names(ds[12]), v2=names(ds[,23:ncol(ds)]))
tmp_ds <- ds[,c(12, 23:ncol(ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(kwFunc, vars[,2], vars[,1], MoreArgs=list(data=tmp_ds), SIMPLIFY=FALSE))
kruskal_wallis_result <- corrs[2:3]
write.csv(kruskal_wallis_result, "metabone-R/results/gender/females/kw/FnStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# kruskal_wallis_result$p.value <- p.adjust(kruskal_wallis_result$p.value, method = "bonferroni", n = length(kruskal_wallis_result$p.value))
# write.csv(kruskal_wallis_result, "metabone-R/results/gender/females/kw/FnStatusvsMetabolitesAdjusted.csv", row.names = TRUE)

## Apply Kruskall-Wallis on TH Status [14]
vars <- data.frame(v1=names(ds[14]), v2=names(ds[,23:ncol(ds)]))
tmp_ds <- ds[,c(14, 23:ncol(ds))]
tmp_ds <- sapply(tmp_ds, as.numeric)
# Apply corrFunc to all rows of vars
corrs <- do.call(rbind, mapply(kwFunc, vars[,2], vars[,1], MoreArgs=list(data=tmp_ds), SIMPLIFY=FALSE))
kruskal_wallis_result <- corrs[2:3]
write.csv(kruskal_wallis_result, "metabone-R/results/gender/females/kw/ThStatusvsMetabolites.csv", row.names = TRUE)
# Apply bonferroni correction ???
# kruskal_wallis_result$p.value <- p.adjust(kruskal_wallis_result$p.value, method = "bonferroni", n = length(kruskal_wallis_result$p.value))
# write.csv(kruskal_wallis_result, "metabone-R/results/gender/females/kw/ThStatusvsMetabolitesAdjusted.csv", row.names = TRUE)


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

## Status Fractures vs Acetate ############################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = Acetate, na.rm = TRUE, fill = Status_Fractures)) +
#  geom_boxplot( , na.rm = TRUE) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Acetate, c(0, 0.9) * 1.1, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "Acetate") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs Acetate", x = "Status_Fractures", y = "Acetate")
p

kruskal.test(ds$Acetate, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(Acetate ~ Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs Lactate ############################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = Lactate, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Lactate, c(0, 0.9) * 1.05, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "Lactate") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs Lactate", x = "Status_Fractures", y = "Lactate")
p

kruskal.test(ds$Lactate, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(Lactate ~ Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs GlycA ##############################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = GlycA, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$GlycA, c(0, 0.9) * 1.1, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "GlycA") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs GlycA", x = "Status_Fractures", y = "GlycA")
p

kruskal.test(ds$GlycA, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(GlycA ~ Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs His ################################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = His, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$His, c(0, 0.9) * 1.1, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "His") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs His", x = "Status_Fractures", y = "His")
p

kruskal.test(ds$His, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(His ~ Status_Fractures, data = ds, exact = FALSE)

## Status Fractures vs Gly ################################################################################
p <- ggplot(ds, aes(x = Status_Fractures, y = Gly, na.rm = TRUE, fill = Status_Fractures)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Gly, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "Status_Fractures", y = "Gly") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_Fractures vs Gly", x = "Status_Fractures", y = "Gly")
p

kruskal.test(ds$Gly, ds$Status_Fractures, rm.na=TRUE)
wilcox.test(Gly ~ Status_Fractures, data = ds, exact = FALSE)


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
## Status BMD vs Val ######################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = Val, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Val, c(0, 0.9) * 1.1, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "Val") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs Val", x = "Status_BMD", y = "Val")
p

by(ds$Val, ds$Status_BMD, shapiro.test)
kruskal.test(ds$Val, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$Val, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## Status BMD vs Citrate ##################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = Citrate, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Citrate, c(0, 0.9) * 1.05, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "Citrate") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs Citrate", x = "Status_BMD", y = "Citrate")
p

by(ds$Citrate, ds$Status_BMD, shapiro.test)
kruskal.test(ds$Citrate, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$Citrate, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## Status BMD vs Leu ######################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = Leu, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Leu, c(0, 0.9) * 1.11, na.rm = TRUE)) +
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

## Status BMD vs TotalBCAA ################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = TotalBCAA, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$TotalBCAA, c(0, 0.9) * 1.10, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "TotalBCAA") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs TotalBCAA", x = "Status_BMD", y = "TotalBCAA")
p

by(ds$TotalBCAA, ds$Status_BMD, shapiro.test)
kruskal.test(ds$TotalBCAA, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$TotalBCAA, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## Status BMD vs His ######################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = His, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$His, c(0, 0.9) * 1.11, na.rm = TRUE)) +
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

## Status BMD vs Ile ######################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = Ile, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Ile, c(0, 0.9) * 1.05, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "Ile") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs Ile", x = "Status_BMD", y = "Ile")
p

by(ds$Ile, ds$Status_BMD, shapiro.test)
kruskal.test(ds$Ile, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$Ile, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## Status BMD vs Phe ######################################################################################
p <- ggplot(ds, aes(x = Status_BMD, y = Phe, na.rm = TRUE, fill = Status_BMD)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Phe, c(0, 0.9) * 1.06, na.rm = TRUE)) +
  labs(x = "Status_BMD", y = "Phe") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "Status_BMD vs Phe", x = "Status_BMD", y = "Phe")
p

by(ds$Phe, ds$Status_BMD, shapiro.test)
kruskal.test(ds$Phe, ds$Status_BMD, rm.na=TRUE)
pairwise.wilcox.test(ds$Phe, ds$Status_BMD, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")


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
## LS status vs LA #######################################################################################
p <- ggplot(ds, aes(x = LS_status, y = LA, na.rm = TRUE, fill = LS_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$LA, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "LS_status", y = "LA") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "LS_status vs LA", x = "LS_status", y = "LA")
p

by(ds$LA, ds$LS_status, shapiro.test)
kruskal.test(ds$LA, ds$LS_status, rm.na=TRUE)
pairwise.wilcox.test(ds$LA, ds$LS_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## LS status vs Citrate ###################################################################################
p <- ggplot(ds, aes(x = LS_status, y = Citrate, na.rm = TRUE, fill = LS_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Citrate, c(0, 0.9) * 1.08, na.rm = TRUE)) +
  labs(x = "LS_status", y = "Citrate") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "LS_status vs Citrate", x = "LS_status", y = "Citrate")
p

by(ds$Citrate, ds$LS_status, shapiro.test)
kruskal.test(ds$Citrate, ds$LS_status, rm.na=TRUE)
pairwise.wilcox.test(ds$Citrate, ds$LS_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")


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
## TH status vs VLDLsize ##################################################################################
p <- ggplot(ds, aes(x = TH_status, y = VLDLsize, na.rm = TRUE, fill = TH_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$VLDLsize, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "TH_status", y = "VLDLsize") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "TH_status vs VLDLsize", x = "TH_status", y = "VLDLsize")
p

by(ds$VLDLsize, ds$TH_status, shapiro.test)
kruskal.test(ds$VLDLsize, ds$TH_status, rm.na=TRUE)
pairwise.wilcox.test(ds$VLDLsize, ds$TH_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")

## TH status vs Phosphatidylc #############################################################################
p <- ggplot(ds, aes(x = TH_status, y = Phosphatidylc, na.rm = TRUE, fill = TH_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$Phosphatidylc, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "TH_status", y = "Phosphatidylc") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "TH_status vs Phosphatidylc", x = "TH_status", y = "Phosphatidylc")
p

by(ds$Phosphatidylc, ds$TH_status, shapiro.test)
kruskal.test(ds$Phosphatidylc, ds$H_status, rm.na=TRUE)
pairwise.wilcox.test(ds$Phosphatidylc, ds$TH_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")


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
## FN status vs VLDLsize ###############################################################################
p <- ggplot(ds, aes(x = FN_status, y = VLDLsize, na.rm = TRUE, fill = FN_status)) +
  geom_boxplot( , na.rm = TRUE, outlier.shape = NA) + coord_cartesian(ylim = quantile(ds$VLDLsize, c(0, 0.9) * 1.11, na.rm = TRUE)) +
  labs(x = "FN_status", y = "VLDLsize") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line("darkgrey")) +
  labs(title = "FN_status vs VLDLsize", x = "FN_status", y = "VLDLsize")
p

by(ds$VLDLsize, ds$FN_status, shapiro.test)
kruskal.test(ds$VLDLsize, ds$FN_status, rm.na=TRUE)
pairwise.wilcox.test(ds$VLDLsize, ds$FN_status, p.adjust.method = "bonferroni", paired = FALSE, alternative= "two.sided")
