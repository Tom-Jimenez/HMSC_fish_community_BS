################################################################################
##### Script 2 : Fit models ----------------------------------------------------
##### Tom Jimenez, 16-11-2026 --------------------------------------------------
################################################################################

remove(list=ls())
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
species <- c(colnames(data_species))
nb_samples_per_species <- apply(data_species, 2, function(x) sum(x > 0))
selected_species <- colnames(data_species)[nb_samples_per_species >= 130]
selected_species

# Modification traits
library(dplyr)
library(stringr)

str(fish_traits)

# Feeding mode : 
# - merge generalist + piscivorous
# - remove herbivorous
fish_traits <- fish_traits %>%
  mutate(feeding.mode = as.character(feeding.mode), 
         feeding.mode = ifelse(feeding.mode %in% c("generalist", "piscivorous"),
                               "generalist/piscivorous", feeding.mode)) %>%
  filter(feeding.mode != "herbivorous") %>%
  mutate(feeding.mode = as.factor(feeding.mode))

# Habitat : 
# - remove bathydemersal and bathypelagic
# - rename reef-associated to demersal
fish_traits <- fish_traits %>%
  filter(!habitat %in% c("bathydemersal", "bathypelagic")) %>%
  mutate(habitat = as.character(habitat),
         habitat = ifelse(habitat == "reef-associated", "demersal", habitat),
         habitat = as.factor(habitat))

# Body shape : 
# - remove "compressiform", "eel-like" and "short and/or deep" species
fish_traits <- fish_traits %>%
  mutate(body.shape = as.character(body.shape)) %>%
  filter(!body.shape %in% c("compressiform", "short and/or deep")) %>%
  mutate(body.shape = as.factor(body.shape))

summary(fish_traits$feeding.mode)
summary(fish_traits$habitat)
summary(fish_traits$body.shape)
summary(fish_traits$spawning.type)
summary(fish_traits$fin.shape)

# Synchronize selected_species with fish_traits

selected_species <- selected_species[selected_species %in% fish_traits$taxon]
selected_species

# We selected environmental covariates.

XData = data.frame(cell = data$cell, 
                   depth=data$depth, 
                   T_floor = data$T_floor,
                   Sal_floor = data$Sal_floor,
                   O2_floor = data$O2_floor,
                   season = as.factor(data$season),
                   chl = data$chl)

# We truncate the species data to presence-absence.

data_species <- data[, colnames(data) %in% selected_species]
Y = as.matrix(data_species)>0
Y = apply(Y,MARGIN = 2,FUN = as.numeric)
dim(Y)

# Species traits

TrData = data.frame(taxon = fish_traits$taxon, 
                    habitat=fish_traits$habitat,
                    body.shape=fish_traits$body.shape,
                    #fin.shape=fish_traits$fin.shape,
                    #spawning.type=fish_traits$spawning.type,
                    feeding.mode=fish_traits$feeding.mode,
                    log_age.maturity = fish_traits$log_age.maturity,
                    #log_length.max = fish_traits$log_length.max,
                    log_offspring.size = fish_traits$log_offspring.size,
                    log_fecundity = fish_traits$log_fecundity,
                    log_growth = fish_traits$log_growth
                    )

# Only the species selected are kept in TrData

library(dplyr)
TrData <- TrData %>%
  filter(taxon %in% selected_species)

summary(TrData$feeding.mode)
summary(TrData$habitat)
summary(TrData$body.shape)
#save(TrData, file = "./data/TrData.RData")

# Check correlations between traits

library(ggcorrplot)
library(ggplot2)
library(reshape2)
library(vcd)
library(RColorBrewer)

names(TrData)

# Numerical traits

traits_num <- c("log_age.maturity", "log_offspring.size", "log_fecundity", "log_growth")

N <- TrData[, traits_num]

cor_x <- cor(N, method = "pearson", use = "pairwise.complete.obs")
cor_x_round <- round(cor_x, 2)

ggcorrplot(
  cor_x_round,
  hc.order = FALSE,
  type = "lower",
  lab = TRUE,
  lab_size = 4,
  show.diag = TRUE,
  outline.color = "white",
  colors = c("#2166AC", "white", "#E31B23")
) +
  scale_fill_gradient2(
    name = "Pearson r",
    low = "#2166AC",
    mid = "white",
    high = "#E31B23",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Pearson pairwise correlation coefficients \n for the numerical traits included in the model.",
    subtitle = "Red (blue) colour indicates positive (negative) coefficients.",
    x = NULL,
    y = NULL
  )

# Categorical traits

traits_cat <- c("habitat", "body.shape", "feeding.mode")
# "fin.shape", "spawning.type"

C <- TrData[, traits_cat]

# Function: Cramer’s V between two factors
cramer_v <- function(x, y) {
  tbl <- table(droplevels(x), droplevels(y))
  if (nrow(tbl) < 2 || ncol(tbl) < 2) return(0)
  assocstats(tbl)$cramer
}

# Cramer’s V-matrix (all pairs of variables)
vars <- names(C)
n    <- length(vars)

mat <- matrix(NA, nrow = n, ncol = n, dimnames = list(vars, vars))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    mat[i, j] <- if (i == j) 1 else cramer_v(C[[i]], C[[j]])
  }
}

mat_df <- melt(mat, varnames = c("Var1", "Var2"), value.name = "CramerV")

mat_df <- mat_df %>%
  mutate(across(c(Var1, Var2), ~ factor(gsub("\\.", " ", as.character(.)))))

vars <- unique(c(as.character(mat_df$Var1), as.character(mat_df$Var2)))

mat_cramer <- matrix(NA_real_, nrow = length(vars), ncol = length(vars),
                     dimnames = list(vars, vars))

for (i in seq_len(nrow(mat_df))) {
  v1 <- as.character(mat_df$Var1[i])
  v2 <- as.character(mat_df$Var2[i])
  val <- mat_df$CramerV[i]
  
  mat_cramer[v1, v2] <- val
  mat_cramer[v2, v1] <- val
}

diag(mat_cramer) <- 1

mat_cramer

ggcorrplot(
  mat_cramer,
  hc.order = FALSE,
  type = "lower",
  lab = TRUE,
  lab_size = 4,
  show.diag = TRUE,
  outline.color = "white",
  colors = c("white", "#E31B23")
) +
  scale_fill_gradient2(
    name = "Cramér's V",
    low = "white",
    high = "#E31B23",
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title    = "Associations between qualitative traits \n measured using Cramér's V",
    subtitle = "Based on the chi-squared test",
    x = NULL,
    y = NULL
  )

# Extracting taxonomic data on selected species

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
tree_plot$edge.length <- rep(1,nrow(tree_plot$edge))
#save(tree_plot, file = "./data/tree_plot.RData")

# Check species order
stopifnot(setequal(tree_plot$tip.label, colnames(Y)))

# Random effects

studyDesign = data.frame(cell = XData$cell, Year = as.factor(data$Year)) 
studyDesign <- studyDesign %>% mutate(across(everything(), as.factor))

# Spatial random effect

# Lon_c / Lat_c are the coordinates of the cell center 
# they are therefore identical for all observations sharing the same cell. 
# One can extract one row per cell directly.

xy_cell_df <- data %>%
  distinct(cell, Lon_c, Lat_c) %>%
  arrange(cell)

# Check: a single pair of coordinates per cell 
# if this is not the case, it indicates a data inconsistency that must be corrected before proceeding
stopifnot(nrow(xy_cell_df) == n_distinct(data$cell))

xy_cell <- as.matrix(xy_cell_df[, c("Lon_c", "Lat_c")])
rownames(xy_cell) <- xy_cell_df$cell
colnames(xy_cell) <- c("x-coordinate", "y-coordinate")

# Check: each level of studyDesign$cell must have an entry in xy_cell
stopifnot(all(levels(studyDesign$cell) %in% rownames(xy_cell)))

rL = HmscRandomLevel(sData = xy_cell, longlat = TRUE,
                     sMethod = "NNGP", nNeighbours = 10)
rL <- setPriors(rL, nfMin = 1, nfMax = 6)
rL

# Temporal random effect

# Temporal autocorrelation structure: years as a 1D coordinate,
# which results in covariance that decreases with temporal distance.

year_levels <- sort(unique(as.numeric(as.character(data$Year))))
yr_coords <- matrix(year_levels, ncol = 1,
                    dimnames = list(as.character(year_levels), "Year"))

rLt = HmscRandomLevel(sData = yr_coords, longlat = FALSE)
rLt <- setPriors(rLt, nfMin = 1, nfMax = 6)
rLt

# Formula

XFormula =  ~  poly(T_floor, degree = 2, raw = TRUE) + 
  poly(Sal_floor, degree = 2, raw = TRUE) +
  poly(depth, degree = 2, raw = TRUE) +
  poly(O2_floor, degree = 2, raw = TRUE) +
  poly(chl, degree = 2, raw = TRUE) +
  season

TrFormula = ~ habitat + 
  body.shape + 
  #fin.shape + 
  #spawning.type + 
  feeding.mode +
  log_age.maturity + 
  #log_length.max + 
  log_offspring.size + 
  log_fecundity +
  log_growth

# Reorder TrData so that it exactly matches the order of species in Y
TrData <- TrData[match(colnames(Y), TrData$taxon), ]

# Check
stopifnot(all(!is.na(TrData$taxon)))
stopifnot(all(TrData$taxon == colnames(Y)))

rownames(TrData) <- TrData$taxon

# We define the model and fit it

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

#save(Hmsc, file = "./data/Hmsc_occurrence.RData")
#load(file = "./data/Hmsc_occurrence.RData")

t0 <- Sys.time()
nChains = 4
nParallel = 4
samples = 500
for (thin in c(500))
{
  transient = round(0.5*samples*thin)
  m = sampleMcmc(Hmsc, thin = thin, samples = samples, transient = transient, #initPar = "fixed effects",
                 nChains = nChains, nParallel = nParallel)
  filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin),"_occurrence"))
  save(m,file=filename)
}
t1 <- Sys.time()
t1-t0
