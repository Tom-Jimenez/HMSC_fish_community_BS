################################################################################
##### Script 2 : Fit models ----------------------------------------------------
################################################################################

# remove(list=ls())
getwd()

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

# We load the Hmsc package, and set the random seed, so that all results from this script are reproductable

library(Hmsc)
set.seed(1)

fish_survey <- readRDS("./data/fish_survey_Baltic_sea_v2.rds")
fish_traits <- readRDS("./data/fish_traits.rds")

fish_survey$Year <- as.factor(fish_survey$Year)
fish_survey$season <- as.factor(fish_survey$season)

# Choose species with more than 130 samples (1% of 12953) in the Baltic Sea

data<-fish_survey
print(colnames(data)[10:150])
data_species <- data[, 10:150]
Species <- c(colnames(data_species))
nb_samples_per_species <- apply(data_species, 2, function(x) sum(x > 0))
selected_species <- colnames(data_species)[nb_samples_per_species >= 130]
selected_species

# Modification traits
library(dplyr)
library(stringr)

str(fish_traits)

# Feeding mode : fusionner generalist + piscivorous
fish_traits <- fish_traits %>%
  mutate(feeding.mode = as.character(feeding.mode), 
         feeding.mode = ifelse(feeding.mode %in% c("generalist", "piscivorous"),
                               "generalist/piscivorous",
                               feeding.mode)) %>%
  filter(feeding.mode != "herbivorous") %>%
  mutate(feeding.mode = as.factor(feeding.mode))     

# Habitat : 
# - supprimer bathydemersal et bathypelagic
# - renommer reef-associated en demersal
fish_traits <- fish_traits %>%
  filter(!habitat %in% c("bathydemersal", "bathypelagic")) %>%
  mutate(habitat = as.character(habitat),
         habitat = ifelse(habitat == "reef-associated", "demersal", habitat),
         habitat = as.factor(habitat))

# Body shape : 
# supprimer les espèces "compressiform" ET "short and/or deep"
fish_traits <- fish_traits %>%
  mutate(body.shape = as.character(body.shape)) %>%
  filter(!body.shape %in% c("compressiform", "short and/or deep")) %>%
  mutate(body.shape = as.factor(body.shape))

summary(fish_traits$feeding.mode)
summary(fish_traits$habitat)
summary(fish_traits$body.shape)

# Synchroniser selected_species avec fish_traits

selected_species <- selected_species[selected_species %in% fish_traits$taxon]
selected_species

# Remove species with higher psrf

species_remove <- c(
  "Molva molva", "Lesueurigobius friesii", "Scophthalmus maximus", "Squalus acanthias")

selected_species <- setdiff(selected_species, species_remove)
selected_species

# We selected to include as candidate environmental covariates only habitat type and spring temperature,
# even if also many other environmental variables could be expected to influence the commmunity.

XData = data.frame(cell = data$cell, 
                   depth=data$depth, 
                   T_floor = data$T_floor,
                   Sal_floor = data$Sal_floor,
                   O2_floor = data$O2_floor,
                   #T_surf = data$T_surf,
                   #Sal_surf = data$Sal_surf,
                   #O2_surf = data$O2_surf,
                   season = as.factor(data$season),
                   chl = data$chl)
                   #ph = data$ph)

# We truncate the species data to presence-absence, and thus ignore abundance variation.

data_species <- data[, colnames(data) %in% selected_species]
Y = as.matrix(data_species)>0
Y = apply(Y,MARGIN = 2,FUN = as.numeric)
dim(Y)

# The data contains also the x- and y-coordinates of the routes. We store these as the xy-matrix to be able to fit a spatial model

xy = as.matrix(cbind(data$lon,data$lat))
rownames(xy)=data$cell
colnames(xy)=c("x-coordinate","y-coordinate")

# Species traits

alltraits = fish_traits
TrData = data.frame(taxon = alltraits$taxon, 
                    habitat=alltraits$habitat,
                    body.shape=alltraits$body.shape,
                    #fin.shape=alltraits$fin.shape,
                    #spawning.type=alltraits$spawning.type,
                    feeding.mode=alltraits$feeding.mode,
                    log_age.maturity = alltraits$log_age.maturity,
                    log_length.max = alltraits$log_length.max,
                    log_offspring.size = alltraits$log_offspring.size,
                    log_fecundity = alltraits$log_fecundity,
                    log_growth = alltraits$log_growth
                    )

# only the species selected from the Baltic Sea are kept in TrData

library(dplyr)
TrData <- TrData %>%
  filter(taxon %in% selected_species)

# Extracting taxonomic data on species present in the Baltic Sea

library(ape)
library(phytools)

load("./data/NS_taxonomy.RData")
taxonomy_selected_species <- my_taxonomy[my_taxonomy$Species %in% selected_species, ]

species_names <- taxonomy_selected_species$Species

taxonomy_selected_species$Phylum <- as.factor(taxonomy_selected_species$Phylum)
taxonomy_selected_species$Class <- as.factor(taxonomy_selected_species$Class)
taxonomy_selected_species$Order <- as.factor(taxonomy_selected_species$Order)
taxonomy_selected_species$Family <- as.factor(taxonomy_selected_species$Family)
taxonomy_selected_species$Genus <- as.factor(taxonomy_selected_species$Genus)
taxonomy_selected_species$Species <- as.factor(taxonomy_selected_species$Species)

taxonomy_selected_species$Path <- apply(taxonomy_selected_species, 1, function(row) {
  paste(row["Phylum"], row["Class"], row["Order"], row["Family"], row["Genus"], row["Species"], sep = "/")
})

tree_plot <- ape::as.phylo.formula(~ Phylum/Class/Order/Family/Genus/Species, data = taxonomy_selected_species)
tree_plot <- as.phylo(tree_plot)

tree_plot$edge.length <- rep(1,length(tree_plot$edge))

# Random effects

studyDesign = data.frame(cell = XData$cell, Year = as.factor(data$Year)) 
studyDesign <- studyDesign %>% mutate(across(everything(), as.factor))

rL = HmscRandomLevel(units = data$cell)
rLt = HmscRandomLevel(units = data$Year)
rL$nfMax = 6
rLt$nfMax = 6

# Formula

XFormula =  ~  poly(T_floor, degree = 2, raw = TRUE) + 
  poly(Sal_floor, degree = 2, raw = TRUE) +
  depth + O2_floor + chl + season
  #T_surf + Sal_surf + O2_surf + ph

TrFormula = ~ habitat + 
  body.shape + 
  #fin.shape + 
  #spawning.type + 
  feeding.mode +
  log_age.maturity + 
  log_length.max + 
  log_offspring.size + 
  log_fecundity +
  log_growth

# Define row names from the ‘taxon’ column
rownames(TrData) = colnames(Y)

# We define the model

Hmsc = Hmsc(Y=Y, 
         XData = XData, 
         XFormula = XFormula, 
         TrData = TrData, 
         TrFormula = TrFormula, 
         phyloTree = tree_plot,
         distr = "probit", 
         studyDesign = studyDesign, 
         ranLevels = list(cell=rL, Year=rLt))

Hmsc

t0 <- Sys.time()
nChains = 4
nParallel = 4
samples = 250
for (thin in c(100))
{
  transient = round(0.5*samples*thin)
  m = sampleMcmc(Hmsc, thin = thin, samples = samples, transient = transient,
                 nChains = nChains, initPar = "fixed effects",
                 nParallel = nParallel)
  filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin),"_occurrence"))
  save(m,file=filename)
}
t1 <- Sys.time()
t1-t0
