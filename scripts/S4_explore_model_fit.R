################################################################################
##### Script 4 : Explore model fit ---------------------------------------------
################################################################################

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")
library(Hmsc)
set.seed(1)

nChains = 4
samples = 250
thin = 100
filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin),"_occurrence"))
#filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin),"_biomass"))
load(filename)

# We next evaluate model fit in terms of explanatory power

preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)
MF

nParallel = 4
run.cross.validation = TRUE # start with TRUE when you introduce the script
filename=file.path(model.directory, paste0("CV_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(m$thin)))
if(run.cross.validation){
  partition = createPartition(m, nfolds = 4, column = "cell")
  preds = computePredictedValues(m,partition=partition, nParallel = nParallel)
  MFCV = evaluateModelFit(hM=m, predY=preds)
  save(partition,MFCV,file=filename)
} else {
  load(filename)
}

tmp = c(mean(MF$AUC), mean(MFCV$AUC), mean(MF$TjurR2), mean(MFCV$TjurR2), mean(MF$RMSE), mean(MFCV$RMSE))
names(tmp)=c("AUC","AUC (CV)","TjurR2","TjurR2 (CV)", "RMSE","RMSE (CV)")
tmp
tmp1 = c(min(MF$AUC), min(MFCV$AUC), min(MF$TjurR2), min(MFCV$TjurR2), min(MF$RMSE), min(MFCV$RMSE))
names(tmp1)=c("AUC","AUC (CV)","TjurR2","TjurR2 (CV)", "RMSE","RMSE (CV)")
tmp1
tmp2 = c(max(MF$AUC), max(MFCV$AUC), max(MF$TjurR2), max(MFCV$TjurR2), max(MF$RMSE), max(MFCV$RMSE))
names(tmp2)=c("AUC","AUC (CV)","TjurR2","TjurR2 (CV)", "RMSE","RMSE (CV)")
tmp2

# Biomass model

preds_biomass = computePredictedValues(m)
MF_biomass = evaluateModelFit(hM=m, predY=preds_biomass)
MF_biomass

nParallel = 4
run.cross.validation = TRUE # start with TRUE when you introduce the script
filename=file.path(model.directory, paste0("CV_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(m$thin),"_biomass"))
if(run.cross.validation){
  partition = createPartition(m, nfolds = 4, column = "cell")
  preds_biomass = computePredictedValues(m,partition=partition, nParallel = 2)
  MFCV_biomass = evaluateModelFit(hM=m, predY=preds_biomass)
  save(partition,MFCV_biomass,file=filename)
} else {
  load(filename)
}

tmp = c(mean(MF_biomass$R2), mean(MFCV_biomass$R2), mean(MF_biomass$RMSE), mean(MFCV_biomass$RMSE))
names(tmp)=c("R2","R2 (CV)", "RMSE","RMSE (CV)")
tmp

# While the code above yields the average results over the species, we may also look at the species-specific results.
# For example, let us plot the explanatory and predictive AUC measures with respect to each other.

par(mfrow=c(1,1))
plot(MF$AUC,MFCV$AUC)
abline(0,1)

# The call to abline(0,1) adds the identity line (y=x) to the plot. For points (=species) below the line,
# explanatory power is greater than predictive power.

Prevalence = colMeans(Y)
plot(Prevalence,MF$AUC)

# If we would have fitted multiple models, we could compare them either based on the cross-validation
# (as was done by above) or with WAIC. WAIC can be computed as follows:

WAIC = computeWAIC(m)
WAIC

# Discriminatory power 

library(dplyr)
library(colorspace)
library(ggplot2)

MF="./data/MF.RData"
load(MF)
MFCV=file.path(model.directory, paste0("CV_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(m$thin)))
load(MFCV)

dPower <- data.frame(Species = rep(m$spNames, times = 2),
                     AUC = c(MF$AUC, MFCV$AUC),
                     TjurR2 = c(MF$TjurR2, MFCV$TjurR2),
                     RMSE = c(MF$RMSE, MFCV$RMSE),
                     Power = rep(c("Explanatory power", "Predictive power"), each = length(m$spNames)))

AUCorder <- dPower%>%
  filter(Power == "Explanatory power")%>%
  arrange(desc(AUC))
Tjurorder <- dPower%>%
  filter(Power == "Explanatory power")%>%
  arrange(desc(TjurR2))
RMSEorder <- dPower%>%
  filter(Power == "Explanatory power")%>%
  arrange(desc(RMSE))

dAUC <- dPower%>%
  arrange(factor(Species, levels = AUCorder$Species))%>%
  mutate(orderCol = rep(c(1:78), each = 2))
dTjur <- dPower%>%
  arrange(factor(Species, levels = Tjurorder$Species))%>%
  mutate(orderCol = rep(c(1:78), each = 2))
dRMSE <- dPower%>%
  arrange(factor(Species, levels = RMSEorder$Species))%>%
  mutate(orderCol = rep(c(1:78), each = 2))

modelPower <- ggpubr::ggarrange(
  ggplot(dAUC, aes(x = orderCol, y = AUC, fill = Power)) + geom_bar(width = 0.8, stat="identity", position = "dodge") +
    labs(x = "", fill = "", tag = "(a)") +
    scale_fill_manual(values = c("black", "grey")) +
    coord_cartesian(ylim = c(0.7, 1)) + 
    theme_bw() +
    theme(legend.position = c(0.85, 0.95)),
  ggplot(dTjur, aes(x = orderCol, y = TjurR2, fill = Power)) + geom_bar(width = 0.8, stat="identity", position = "dodge") +
    labs(x = "", fill = "", tag = "(b)") +
    scale_fill_manual(values = c("black", "grey")) +
    theme_bw() +
    theme(legend.position = c(0.85, 0.95)),
  ggplot(dRMSE, aes(x = orderCol, y = RMSE, fill = Power)) + geom_bar(width = 0.8, stat="identity", position = "dodge") +
    labs(x = "Species", fill = "", tag = "(b)") +
    scale_fill_manual(values = c("black", "grey")) +
    theme_bw() +
    theme(legend.position = c(0.85, 0.95)),
  common.legend = TRUE, ncol = 1, nrow = 3) ; modelPower

# Table with mean results
dTable <- dPower%>%
  group_by(Power)%>%
  summarise(AUC = round(mean(AUC), 3),
            TjurR2 = round(mean(TjurR2), 3),
            RMSE = round(mean(RMSE), 3)) ; dTable

# Filter species with AUC less than 0.7
low_auc_species <- dPower[dPower$AUC < 0.7, ]

# Filter by "Explanatory power" and "Predictive power"
low_auc_explanatory <- subset(low_auc_species, Power == "Explanatory power")
low_auc_predictive <- subset(low_auc_species, Power == "Predictive power")

# Display the results
cat("Species with AUC < 0.7 - Explanatory Power:\n")
print(low_auc_explanatory$Species)

cat("\nSpecies with AUC < 0.7 - Predictive Power:\n")
print(low_auc_predictive$Species)
