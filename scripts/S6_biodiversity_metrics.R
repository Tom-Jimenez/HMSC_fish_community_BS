################################################################################
##### Script 7 : Biodiversity metrics and trends -------------------------------
################################################################################

# remove(list=ls())

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

################################################################################
##### Creation EpredY for all years --------------------------------------------
################################################################################

library(tidyverse)
library(viridis)
library(vioplot)
library(abind)
library(RColorBrewer)
library(ape)
library(corrplot)
library(gridExtra)
library(grid)
library(sf)       
library(dplyr)    
library(readr) 
library(dggridR)
library(Hmsc)
set.seed(1)

nChains = 4
samples = 250
thin = 100
filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin),"_occurrence"))
#filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin),"_biomass"))
load(filename)

data <- readRDS("./data/fish_survey_Baltic_sea_v2.rds")

# Calculate the longitude and latitude limits in `data`.
min_lon_data <- min(data$lon, na.rm = TRUE) - 1
max_lon_data <- max(data$lon, na.rm = TRUE) + 1
min_lat_data <- min(data$lat, na.rm = TRUE) - 1
max_lat_data <- max(data$lat, na.rm = TRUE) + 1

year=1993

grid <- readRDS(paste0("./data/grid_Copernicus_Baltic_1993-2021/grid_Copernicus_Baltic_", year, ".rds"))
grid <- grid[
  grid$lon >= min_lon_data & grid$lon <= max_lon_data &
    grid$lat >= min_lat_data & grid$lat <= max_lat_data, 
]
grid$season <- "Winter"

# Remove points ‘too-far’ from observations
too_far <- mgcv::exclude.too.far(grid$lon, grid$lat, 
                                 data$lon, data$lat, 
                                 dist = 0.01)
grid <- grid[!too_far, ]

xy.grid = as.matrix(cbind(grid$lon,grid$lat))
XData.grid = data.frame(Sal_floor = grid$Sal_floor, T_floor = grid$T_floor, O2_floor = grid$O2_floor,
                        Sal_surf = grid$Sal_surf, T_surf = grid$T_surf, O2_surf = grid$O2_surf,
                        depth = grid$depth, chl = grid$chl, ph = grid$ph, 
                        season = grid$season, stringsAsFactors = TRUE)
Gradient = prepareGradient(m, XDataNew = XData.grid, sDataNew = list(cell=xy.grid))

predY = predict(m, Gradient=Gradient, expected = TRUE, nParallel = 4)

EpredY = Reduce("+", predY) / length(predY)

saveRDS(EpredY, paste0("./data/EpredY/EpredY_", year, ".rds"))
#saveRDS(EpredY, paste0("./data/EpredY/EpredY_", year, "_biomass.rds"))

remove(predY, EpredY)
gc()

# Same wth a loop
# List of years
years <- 1993:2021
# Loop to process each year in grid + Pred and Epred
for (year in years) {
  
  # Charger les données de grid pour l'année en cours
  grid <- readRDS(paste0("./data/grid_Copernicus_Baltic_1993-2021/grid_Copernicus_Baltic_", year, ".rds"))
  
  # Filter `grid` according to `data` limits
  grid <- grid[
    grid$lon >= min_lon_data & grid$lon <= max_lon_data &
      grid$lat >= min_lat_data & grid$lat <= max_lat_data, 
  ]
  #grid <- grid %>% slice(seq(1, n(), by = 10))
  
  # Add a season column to the grid
  grid$season <- "Winter"
  
  # Remove points ‘too-far’ from observations
  too_far <- mgcv::exclude.too.far(grid$lon, grid$lat, 
                                   data$lon, data$lat, 
                                   dist = 0.01)
  grid <- grid[!too_far, ]
  
  xy.grid = as.matrix(cbind(grid$lon,grid$lat))
  
  XData.grid = data.frame(Sal_floor = grid$Sal_floor, T_floor = grid$T_floor, O2_floor = grid$O2_floor,
                          Sal_surf = grid$Sal_surf, T_surf = grid$T_surf, O2_surf = grid$O2_surf,
                          depth = grid$depth, chl = grid$chl, ph = grid$ph, 
                          season = grid$season, stringsAsFactors = TRUE)
 
  Gradient = prepareGradient(m, XDataNew = XData.grid, sDataNew = list(cell=xy.grid))
  
  nParallel = 4
  predY = predict(m, Gradient=Gradient, expected = TRUE, nParallel=nParallel)

  EpredY = Reduce("+", predY) / length(predY)
  
  # Save
  saveRDS(EpredY, paste0("./data/EpredY/EpredY_", year, ".rds"))
  #saveRDS(EpredY, paste0("./data/EpredY/EpredY_", year, "_biomass.rds"))
  
  cat("Calculation of predictions for the year", year, "finished.\n")
  remove(predY, EpredY)
  gc()
}

# Combine all the years
# List of years
years <- 1993:2021
# Initialise an empty DataFrame to store all rows
combined_data <- data.frame()

# Loop to load each EpredY file and add it to the DataFrame
for (year in years) {
  EpredY <- readRDS(paste0("./data/EpredY/EpredY_", year, ".rds"))
  #EpredY <- readRDS(paste0("./data/EpredY/EpredY_", year, "_biomass.rds"))
  EpredY <- as.data.frame(EpredY)
  # Add a ‘year’ column to identify the year of each observation
  EpredY$year <- year
  # Add lines from the current file to the combined DataFrame
  combined_data <- rbind(combined_data, EpredY)
  cat("Data for the year", year, "added.\n")
}

grid <- readRDS("./data/grid_Copernicus_Baltic_1993-2021/grid_Copernicus_Baltic_2021.rds")
grid <- grid[
  grid$lon >= min_lon_data & grid$lon <= max_lon_data &
    grid$lat >= min_lat_data & grid$lat <= max_lat_data, ]

# Remove points ‘too-far’ from observations
too_far <- mgcv::exclude.too.far(grid$lon, grid$lat, 
                                 data$lon, data$lat, 
                                 dist = 0.01)
grid <- grid[!too_far, ]

expanded_grid <- grid[rep(1:nrow(grid), times = 29), ]
expanded_grid <- expanded_grid[1:nrow(combined_data), ]

combined_data$lon <- expanded_grid$lon
combined_data$lat <- expanded_grid$lat
combined_data$X <- NULL
library(dplyr)
combined_data <- combined_data %>%
  select(lon, lat, year, everything())

colnames(combined_data)[4:43] <- gsub("\\.", " ", colnames(combined_data)[4:43])
names(combined_data)[names(combined_data) == "lon"] <- "longitude"
names(combined_data)[names(combined_data) == "lat"] <- "latitude"

saveRDS(combined_data, file = paste0("./data/EpredY/Combined_EpredY_1993-2021.rds"))
#saveRDS(combined_data, file = paste0("./data/EpredY/Combined_EpredY_1993-2021_biomass.rds"))
cat("The combined file has been saved.\n")

# Epredy Average
years <- 1993:2021

# Initialize a running sum DataFrame and a count DataFrame
sum_EpredY <- NULL
count_EpredY <- NULL

# Loop to load each EpredY file and update the sum and count
for (year in years) {
  # Load the EpredY file for the current year
  EpredY <- readRDS(paste0("./data/EpredY/EpredY_", year, ".rds"))
  #EpredY <- readRDS(paste0("./data/EpredY/EpredY_", year, "_biomass.rds"))
  EpredY <- as.data.frame(EpredY)
  # If it's the first year, initialize the sum and count DataFrames
  if (is.null(sum_EpredY)) {
    sum_EpredY <- as.matrix(EpredY)
    count_EpredY <- !is.na(EpredY)
  } else {
    # Add the current year to the sum, ignoring NA values
    sum_EpredY <- sum_EpredY + as.matrix(EpredY, na.rm = TRUE)
    count_EpredY <- count_EpredY + !is.na(EpredY)
  }
  cat("Data for year", year, "processed.\n")
}

# Compute the mean by dividing the sum by the count
EpredY_mean <- sum_EpredY / count_EpredY

# Replace NaN values (from division by zero) with NA
EpredY_mean[is.nan(EpredY_mean)] <- NA

# Convert to DataFrame for easier handling
EpredY_mean <- as.data.frame(EpredY_mean)
colnames(EpredY_mean)[1:40] <- gsub("\\.", " ", colnames(EpredY_mean)[1:40])

min_lon_data <- min(data$lon, na.rm = TRUE) - 1
max_lon_data <- max(data$lon, na.rm = TRUE) + 1
min_lat_data <- min(data$lat, na.rm = TRUE) - 1
max_lat_data <- max(data$lat, na.rm = TRUE) + 1

grid <- readRDS("./data/grid_Copernicus_Baltic_1993-2021/grid_Copernicus_Baltic_2021.rds")
grid <- grid[
  grid$lon >= min_lon_data & grid$lon <= max_lon_data &
    grid$lat >= min_lat_data & grid$lat <= max_lat_data, ]

# Remove points ‘too-far’ from observations
too_far <- mgcv::exclude.too.far(grid$lon, grid$lat, 
                                 data$lon, data$lat, 
                                 dist = 0.01)
grid <- grid[!too_far, ]

EpredY_mean$lon <- grid$lon
EpredY_mean$lat <- grid$lat

library(dplyr)
EpredY_mean <- EpredY_mean %>%
  select(lon, lat, everything())

names(EpredY_mean)[names(EpredY_mean) == "lon"] <- "longitude"
names(EpredY_mean)[names(EpredY_mean) == "lat"] <- "latitude"

saveRDS(EpredY_mean, file = paste0("./data/EpredY/EpredY_mean_1993-2021.rds"))
#saveRDS(EpredY_mean, file = paste0("./data/EpredY/EpredY_mean_1993-2021_biomass.rds"))

################################################################################
##### Biodiversity metrics ----------------------------------------------
################################################################################

library(dplyr)
library(purrr)
library(tidyverse)
library(vegan)
library(mFD)
library(RColorBrewer)
library(knitr)
library(magrittr)
theme_set(theme_bw())
library(ggplot2)
library(patchwork)

# Load predictions from model: probability of occurrence and (log) biomass

p.occurrence <- readRDS("./data/EpredY/Combined_EpredY_1993-2021.rds")
biomass <- readRDS("./data/EpredY/Combined_EpredY_1993-2021_biomass.rds")
#p.occurrence <- readRDS("./data/EpredY/EpredY_mean_1993-2021.rds")
#biomass <- readRDS("./data/EpredY/EpredY_mean_1993-2021_biomass.rds")

p.occurrence[1:5, 1:10] %>% kable()
biomass[1:5, 1:10] %>% kable()

# Are all columns in the same order?  

table(names(p.occurrence) == names(biomass)) #TRUE

all(biomass %>% select(year, longitude, latitude) == p.occurrence %>% select(year, longitude, latitude)) # TRUE

# Create vector with species names

species <- names(p.occurrence)[!names(p.occurrence) %in%
                                 c("year", "longitude", "latitude")]
species[1:10]

##### Taxonomic diversity ######################################################
##### Taxonomic richness -------------------------------------------------------

# Set a threshold of e.g., 0.2 to determine a species presence 
# one could try different thresholds as sensitivity test

threshold <- 0.2
realized_occurrence <- p.occurrence %>%
  mutate_at(.vars = species, # Apply transformation only to species
            .funs = function(x) ifelse(x >= threshold, 1, 0))

realized_occurrence[1:5, 1:10] %>% kable()

# Compute richness, i.e., number of species presences (1) on each row

richness <- cbind(realized_occurrence %>% select(-all_of(species)), # Keep time and location variables
                  Richness = rowSums(realized_occurrence %>% select(all_of(species)))) # Sum of species occurrences

range(richness$Richness)

##### Shannon & evenness -------------------------------------------------------

# Since data is in log scale (i.e., log(x + 1)), back-transform data to compute Shannon (i.e., exp(x) - 1)

biomass[, species] <- exp(biomass[, species]) - 1
biomass[, species][biomass[, species]<0] <- 0
biomass[, species] <- biomass[, species] * p.occurrence[, species]

shannon <- diversity(biomass[, species], index = "shannon")
shannon[1:10]

# To obtain evennes, we need the number of species considered to be present when computing Shannon
# N species present 

nSp <- rowSums(ifelse(biomass[, species] > 0, 1, 0))
nSp[1:10]

# Compute evenness

evennes <- shannon/log(nSp)
all(biomass %>% select(year, longitude, latitude) == richness %>% select(year, longitude, latitude)) # TRUE

# Merge with previous biodiversity metrics and visualize

fish.diversity <- cbind(richness,
                        potentialRichness = nSp,
                        Shannon = shannon,
                        Evenness = evennes) %>%
  mutate(year = as.numeric(as.character(year)))

# Note that richness values are much higher if we consider p(occ) > 0 as a presence, than when setting a threshold (e.g., 0.5)

ggplot(fish.diversity, aes(x = Richness, y = potentialRichness)) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4))

##### Functional diversity #####################################################

# The computation of functional diversity is based on the mFD package
# Essentially there are two ways to obtain different functional biodiversity metrics:
# 1 - Using the Hill numbers
# 2 - Using a multidimensional space based on a PCoA of species traits

##### Hill numbers -------------------------------------------------------------

# Load species traits

traits <- readRDS("./data/fish_traits.rds")
traits[1:5, 1:10] %>% kable()
traits <- traits[traits$taxon %in% species, ]

# Explore possible correlations among trais, here we set a threshold of 0.7

traits %>%
  select_if(is.numeric) %>%
  cor() %>%
  as.data.frame() %>%
  mutate_all(.funs = function(x) round(x, 2)) %>%
  mutate_all(.funs = function(x) ifelse(abs(x) < 0.7, NA, x)) %>%
  kable()

# length max and length infinity have correlations > 0.7. Let’s remove length infinity to avoid correlated traits in the analysis

traits$log_length.infinity <- NULL

traits %>%
  select_if(is.numeric) %>%
  cor() %>%
  as.data.frame() %>%
  mutate_all(.funs = function(x) round(x, 2)) %>%
  mutate_all(.funs = function(x) ifelse(abs(x) < 0.7, NA, x)) %>%
  kable()

# length max have correlations > 0.7 with growth and age max. Let’s remove length max to avoid correlated traits in the analysis

traits$log_length.max <- NULL

traits %>%
  select_if(is.numeric) %>%
  cor() %>%
  as.data.frame() %>%
  mutate_all(.funs = function(x) round(x, 2)) %>%
  mutate_all(.funs = function(x) ifelse(abs(x) < 0.7, NA, x)) %>%
  kable()


# Now let’s remove unnecessary information from the trait dataframe and make sure that species names (i.e., taxon) is set as rownames

traits <- traits %>%
  select(-c(family, genus, species)) %>%
  column_to_rownames("taxon")
traits[1:5, 1:10] %>% kable()

# Retrieve trait names and type (i.e., ordinal, categorical continuous, fuzzy) into a table

tr_cat <- tibble(trait_name  = names(traits),
                 trait_type = c(rep("N", times = 5), # Nominal traits (factor variable)
                                rep("Q", times = 7))) # Quantitative traits (numeric values)
head(tr_cat) %>% kable()

# Compute species distances

sp_dist <- funct.dist(sp_tr = traits,
                      tr_cat = tr_cat, 
                      # Since there are categorical traits, we compute gower distance
                      metric = "gower", 
                      # If only continuous traits, one can use euclidean distance instead
                      # metric = "euclidean",
                      scale_euclid = "scale_center")

as.matrix(sp_dist)[1:5, 1:10] %>% kable()

# Now one can use these functional distances to compute a series of functional diversity metrics (Hill numbers)

realized_occurrence[1:5, 1:10] %>% kable()

rownames(realized_occurrence) <- NULL
pa.assamblage <- realized_occurrence[, species] %>%
  mutate(ID = paste("community", 1:nrow(realized_occurrence), sep = "_")) %>%
  column_to_rownames("ID") %>%
  as.matrix()

pa.assamblage[1:5, 1:10] %>% kable()

# Compute Hill numbers

hill_test <- alpha.fd.hill(asb_sp_w = pa.assamblage[1:1000,],
                           sp_dist = sp_dist,
                           tau = "mean",
                           q = 0) # Richness-like
head(hill_test$asb_FD_Hill) %>% kable()

hill0 <- alpha.fd.hill(asb_sp_w = pa.assamblage,
                       sp_dist = sp_dist,
                       tau = "mean",
                       q = 0) # Richness-like

rownames(biomass) <- NULL
biomass.assamblage <- biomass[, species] %>%
  mutate(ID = paste("community", 1:nrow(realized_occurrence), sep = "_")) %>%
  column_to_rownames("ID") %>%
  as.matrix()

biomass.assamblage[1:5, 1:10] %>% kable()

hill12 <- alpha.fd.hill(asb_sp_w = biomass.assamblage,
                        sp_dist = sp_dist,
                        tau = "mean",
                        q = c(1, 2)) # Shannon-like & # Simpson-like

# Make sure that the data used when computing Hill numbers has same order than other biodiversity metrics, and merge

table(fish.diversity %>% select(year, longitude, latitude) == realized_occurrence %>% select(year, longitude, latitude))

fish.diversity <- cbind(fish.diversity,
                        funHill0 = hill0$asb_FD_Hill[, 1],
                        funHill1 = hill12$asb_FD_Hill[, 1],
                        funHill2 = hill12$asb_FD_Hill[, 2])

##### Multidimensional functional diversity ------------------------------------

# Functional diversity metrics based on a PCoA (multidimensional)
# The first step is to identify the PCoA settings that best depict the community: 
# - Maximum number of axis (maxdim_pcoa) 
# - Method used to weight the differences between species pairwise distance (deviation_weighting): squared (“mad” in output) or absolute (“rmsd”) 
# - Scaling (fdist_scaling)

quality_fspaces <- quality.fspaces(sp_dist = sp_dist,
                                   maxdim_pcoa = 15,
                                   fdendro = "average",
                                   deviation_weighting = c("absolute", "squared"),
                                   fdist_scaling = c(TRUE, FALSE))

head(quality_fspaces$quality_fspaces) %>% kable()

# Identify the functional space associated with minimal quality metric:
# number of dimensions where the quality is best (smaller value)

apply(quality_fspaces$quality_fspaces, 2, which.min) 

# The actual value

apply(quality_fspaces$quality_fspaces, 2, min) 

# mad_scaled with 5 dimensions has the smallest value => this are the best settings to represent the community from a functional perspective
# We can also visualize the quality of the selected functional spaces

quality.fspaces.plot(fspaces_quality = quality_fspaces,
                     quality_metric = "mad_scaled",
                     fspaces_plot = c("tree_average", "pcoa_6d", "pcoa_7d", "pcoa_8d"),
                     # fspaces_plot = c("tree_average", "pcoa_2d", "pcoa_3d", "pcoa_4d", "pcoa_5d", "pcoa_6d", "pcoa_7d", "pcoa_8d", "pcoa_9d", "pcoa_10d"),
                     gradient_deviation = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
                     gradient_deviation_quality = c(low = "yellow", high = "red"),
                     x_lab = "Trait-based distance")

# For the 2D space, on the top row there are a lot of points below the 1:1 lines, meaning that distances are overestimated in this multidimensional space. 

quality_fspaces$"quality_fspaces" %>%
  tibble::as_tibble(rownames = "Funct.space") %>%
  tidyr::pivot_longer(cols =! Funct.space, names_to = "quality_metric", values_to = "Quality") %>%
  ggplot(aes(x = Funct.space, y = Quality, color = quality_metric, shape = quality_metric)) +
  geom_point() + lims(y = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# As described previously, the PCoA with mad_scaled and 8 dimensions is the setup that best describes the community from a functional perspective. Let’s proceed with that

mad_scaled_fspaces <- quality.fspaces(sp_dist = sp_dist,
                                      maxdim_pcoa = 8,
                                      fdendro = "average",
                                      deviation_weighting = "absolute",
                                      fdist_scaling = TRUE)

# Even if the PCoA has 8 dimensions, might be a good idea to use a smaller number of PCs to compute the functional metrics, since the more PCs, the longer the computation takes. 
# We can thus identify the number of PCs that explains e.g., 75% of the variance.

variance_threshold <- 0.75
eigenvalues <- mad_scaled_fspaces$details_fspaces$pc_eigenvalues$Eigenvalues

ggplot(data.frame(Proportion = cumsum(eigenvalues)/sum(eigenvalues),
                  Number_axis = 1:length(eigenvalues)),
       aes(x = Number_axis, y = Proportion)) +
  geom_abline(intercept = variance_threshold, slope = 0, color = "red") +
  geom_line() +
  geom_point()

# To explain >75% of variance one should use 5 PCoA axis in the calculations

sp_faxes_coord <- mad_scaled_fspaces$"details_fspaces"$"sp_pc_coord"

tr_faxes_Q <- traits.faxes.cor(sp_tr = traits[, tr_cat$trait_name[tr_cat$trait_type == "Q"]], # Only QUANTITATIVE traits
                               sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3")], # Display only relationships with the first 3 axis (PCs)
                               plot = TRUE) # return plot
tr_faxes_Q

tr_faxes_N <- traits.faxes.cor(sp_tr = traits[, tr_cat$trait_name[tr_cat$trait_type == "N"]], # Only CATECORICAL traits
                               sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3")],
                               plot = TRUE)
tr_faxes_N

# Finally we can compute the functional diversity metrics based on the PCoA. A reminder that even though the PCoA has 8 PCs, we only use 5 to slow down computation

biomass.assamblage[1:5, 1:10] %>% kable()

# Test with few observations first

test_fundiv <- alpha.fd.multidim(
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5")],
  asb_sp_w = biomass.assamblage[1:100, ],
  # Define which functional diversity metrics are to be computed
  ind_vect = c( "fric", "feve", "fdis", "fdiv"),
  scaling = TRUE,
  check_input = TRUE,
  details_returned = TRUE)

head(test_fundiv$functional_diversity_indices) %>% kable()

# Running the whole set of communities (i.e., rows) 

alpha_fd_indices <- alpha.fd.multidim(
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5")],
  asb_sp_w = biomass.assamblage,
  # Define which functional diversity metrics are to be computed
  # Note that more indices will slow down the process
  ind_vect = c( "fric", "feve", "fdis", "fdiv"),
  scaling = TRUE,
  check_input = TRUE,
  details_returned = TRUE)

head(alpha_fd_indices$functional_diversity_indices) %>% kable()

# Again, for computing functional richness it might be more adequate to use only presence/absence data

pa.assamblage[1:5, 1:10] %>% kable()

# It also takes quite long, so here I provide the already computed values

f_ric <- alpha.fd.multidim(
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5")],
  asb_sp_w = pa.assamblage,
  # Define which functional diversity metrics are to be computed
  ind_vect = c( "fric"),
  scaling = TRUE,
  check_input = TRUE,
  details_returned = TRUE)

head(f_ric$functional_diversity_indices) %>% kable()

# We can finally merge the multivariate functional diversity metrics with all other

fish.diversity <- cbind(fish.diversity,
                        fric_pa = f_ric$functional_diversity_indices$fric,
                        alpha_fd_indices$functional_diversity_indices)

fish.diversity[1:5, 1:10] %>% kable()

# Another way of calculating indices by subset

library(tidyverse)
library(dplyr)

# Calculation of hill0 by subsets 

n_rows <- nrow(pa.assamblage)  
chunk_size <- 100000          
output_dir <- "./data/Biodiversity metrics/output_hill0_rds"  
dir.create(output_dir, showWarnings = FALSE)  

for (start_row in seq(1, n_rows, by = chunk_size)) {
  end_row <- min(start_row + chunk_size - 1, n_rows) 
  message(paste("Processing of lines", start_row, "to", end_row))
  hill_result <- alpha.fd.hill(
    asb_sp_w = pa.assamblage[start_row:end_row, ],
    sp_dist = sp_dist,
    tau = "mean",
    q = 0)
  output_file <- file.path(output_dir, paste0("hill0_", start_row, "_", end_row, ".rds"))
  saveRDS(hill_result, file = output_file)
}

# Combine results in a single file

hill_combined <- map_dfr(
  list.files(output_dir, full.names = TRUE, pattern = "\\.rds$"),
  ~ {
    data <- readRDS(.x)
    tibble(
      community = rownames(data$asb_FD_Hill),
      FD_q0 = data$asb_FD_Hill[, 1],
      tau_dist = data$tau_dist,
      asb_totw = data$details$asb_totw
    )
  }
)

# Reorganisation of hill0 lines according to the ‘community’ column

hill_combined <- hill_combined %>%
  mutate(community_num = as.numeric(gsub("community_", "", community))) %>%
  arrange(community_num) %>%
  select(-community_num) 

saveRDS(hill_combined, "./data/Biodiversity metrics/Hill0_combined.rds")

# Calculation of hill12 by subsets 

n_rows <- nrow(biomass.assamblage)
chunk_size <- 100000      
output_dir <- "./data/Biodiversity metrics/output_hill12_rds" 
dir.create(output_dir, showWarnings = FALSE) 

for (start_row in seq(1, n_rows, by = chunk_size)) {
  end_row <- min(start_row + chunk_size - 1, n_rows)
  message(paste("Processing of lines", start_row, "to", end_row))
  hill_result <- alpha.fd.hill(
    asb_sp_w = biomass.assamblage[start_row:end_row, ],
    sp_dist = sp_dist,
    tau = "mean",
    q = c(1, 2))
  output_file <- file.path(output_dir, paste0("hill12_", start_row, "_", end_row, ".rds"))
  saveRDS(hill_result, file = output_file)
}

hill_combined <- map_dfr(
  list.files(output_dir, full.names = TRUE, pattern = "\\.rds$"),
  ~ {
    data <- readRDS(.x)
    tibble(
      community = rownames(data$asb_FD_Hill),
      FD_q1 = data$asb_FD_Hill[, 1],
      FD_q2 = data$asb_FD_Hill[, 2],
      tau_dist = data$tau_dist,
      asb_totw = data$details$asb_totw
    )
  }
)

hill_combined <- hill_combined %>%
  mutate(community_num = as.numeric(gsub("community_", "", community))) %>%
  arrange(community_num) %>%
  select(-community_num)

saveRDS(hill_combined, "./data/Biodiversity metrics/Hill12_combined.rds")

hill0 <- read_rds("./data/Biodiversity metrics/Hill0_combined.rds")
hill12 <- read_rds("./data/Biodiversity metrics/Hill12_combined.rds")
fish.diversity <- cbind(fish.diversity,
                        funHill0 = hill0$FD_q0,
                        funHill1 = hill12$FD_q1,
                        funHill2 = hill12$FD_q2)

# Calculation of alpha_fd_indices by subsets 

n_rows <- nrow(biomass.assamblage)
chunk_size <- 10000      
output_dir <- "./data/Biodiversity metrics/output_alpha_fd_indices_rds" 
dir.create(output_dir, showWarnings = FALSE) 

for (start_row in seq(1, n_rows, by = chunk_size)) {
  end_row <- min(start_row + chunk_size - 1, n_rows)
  message(paste("Processing of lines", start_row, "to", end_row))
  fundiv_result <- alpha.fd.multidim(
    sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5")],
    asb_sp_w = biomass.assamblage[start_row:end_row, ],
    ind_vect = c( "fric", "feve", "fdis", "fdiv"),
    scaling = TRUE,
    check_input = F,
    details_returned = F)
  output_file <- file.path(output_dir, paste0("alpha_fd_indices_", start_row, "_", end_row, ".rds"))
  saveRDS(fundiv_result, file = output_file)
}

fundiv_combined <- map_dfr(
  list.files(output_dir, full.names = TRUE, pattern = "\\.rds$"),
  ~ {
    data <- readRDS(.x)
    tibble(
      community = rownames(data$functional_diversity_indices),
      data$functional_diversity_indices
    )
  }
)

fundiv_combined <- fundiv_combined %>%
  mutate(community_num = as.numeric(gsub("community_", "", community))) %>%
  arrange(community_num) %>%
  select(-community_num)

saveRDS(fundiv_combined, "./data/Biodiversity metrics/alpha_fd_indices_combined.rds")

# Calculation of f_ric by subsets

n_rows <- nrow(pa.assamblage)
chunk_size <- 10000      
output_dir <- "./data/Biodiversity metrics/output_f_ric_rds" 
dir.create(output_dir, showWarnings = FALSE) 

for (start_row in seq(1, n_rows, by = chunk_size)) {
  end_row <- min(start_row + chunk_size - 1, n_rows)
  message(paste("Processing of lines", start_row, "to", end_row))
  f_ric_result <- alpha.fd.multidim(
    sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5")],
    asb_sp_w = pa.assamblage[start_row:end_row, ],
    ind_vect = c( "fric"),
    scaling = TRUE,
    check_input = F,
    details_returned = F)
  output_file <- file.path(output_dir, paste0("f_ric_", start_row, "_", end_row, ".rds"))
  saveRDS(f_ric_result, file = output_file)
}

f_ric_combined <- map_dfr(
  list.files(output_dir, full.names = TRUE, pattern = "\\.rds$"),
  ~ {
    data <- readRDS(.x)
    tibble(
      community = rownames(data$functional_diversity_indices),
      fric_pa = data$functional_diversity_indices$fric
    )
  }
)

f_ric_combined <- f_ric_combined %>%
  mutate(community_num = as.numeric(gsub("community_", "", community))) %>%
  arrange(community_num) %>%
  select(-community_num)

saveRDS(f_ric_combined, "./data/Biodiversity metrics/f_ric_combined.rds")

alpha_fd_indices <- read_rds("./data/Biodiversity metrics/alpha_fd_indices_combined.rds")
fish.diversity <- cbind(fish.diversity,
                        fric = alpha_fd_indices$fric,
                        feve = alpha_fd_indices$feve,
                        fdis = alpha_fd_indices$fdis,
                        fdiv = alpha_fd_indices$fdiv,
                        fide_PC1 = alpha_fd_indices$fide_PC1,
                        fide_PC2 = alpha_fd_indices$fide_PC2,
                        fide_PC3 = alpha_fd_indices$fide_PC3,
                        fide_PC4 = alpha_fd_indices$fide_PC4,
                        fide_PC5 = alpha_fd_indices$fide_PC5)

f_ric <- read_rds("./data/Biodiversity metrics/f_ric_combined.rds")
fish.diversity <- cbind(fish.diversity,
                        fric_pa = f_ric$fric_pa
)
#fish.diversity$community <- NULL

# saveRDS(fish.diversity, file = "./data/Biodiversity metrics/fish.diversity.rds")
# fish.diversity <- read_rds("./data/Biodiversity metrics/fish.diversity.rds")
# saveRDS(biomass.assamblage, file = "./data/Biodiversity metrics/biomass.assamblage.rds")
# biomass.assamblage <- read_rds("./data/Biodiversity metrics/biomass.assamblage.rds")
# saveRDS(pa.assamblage, file = "./data/Biodiversity metrics/pa.assamblage.rds")
# pa.assamblage <- read_rds("./data/Biodiversity metrics/pa.assamblage.rds")

# PCA ----
library(FactoMineR)
library(factoextra)

fish.diversity_mean <- read_rds("./data/Biodiversity metrics/fish.diversity_mean.rds")
selected_data <- fish.diversity_mean[, c("Richness", "Shannon", "funHill0", "fric_pa", "fdiv")]
colnames(selected_data) <- c(
  "Richness",
  "Shannon",
  "Hill 0",
  "Functional richness",
  "Functional divergence"
)

selected_data <- fish.diversity_mean[, c("Trend Richness", "Trend Shannon", "Trend funHill0", "Trend fric_pa", "Trend fdiv")]
colnames(selected_data) <- c(
  "Richness",
  "Shannon",
  "Hill 0",
  "Functional richness",
  "Functional divergence"
)

res_pca <- PCA(selected_data, graph = FALSE)

p <- fviz_pca_var(res_pca,
                  col.var = "contrib", 
                  gradient.cols = c('blue', "red", "orange"),
                  repel = TRUE, 
                  axes = c(1, 2)) 
p <- p + theme(
  legend.position = c(0.15, 0.25),
  legend.background = element_rect(fill = "white", color = "white")
) +
  annotate("text", x = -Inf, y = Inf, label = "L", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")  
p

# Plots

fish.diversity <- read_rds("./data/Biodiversity metrics/fish.diversity.rds")
fish.diversity_mean <- read_rds("./data/Biodiversity metrics/fish.diversity_mean.rds")

plotFun <- function(data, variable, years = NULL, size = 0.5){
  if(!is.null(years)){
    data <- data %>%
      filter(year %in% years)
  }
  p <- ggplot(data, aes(x = longitude, y = latitude, color = .data[[variable]])) +
    geom_point(size = size) + 
    borders(fill = "grey", colour = "black", size = 0.25) + 
    coord_quickmap(xlim = range(data$longitude), ylim = range(data$latitude)) +
    scale_color_gradientn(colours = rev(brewer.pal(11, "RdYlBu"))) +
    annotate("text", x = min(data$longitude), y = max(data$latitude), label = "J", hjust = 0, vjust = 1, size = 5) +
    theme(
      legend.position = c(0.95, 0.15),
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_blank()  
    )
  
  if(!is.null(years)){
    p <- p + facet_wrap(~ year)
  }
  return(p)
}

# Compute linear trend over time
r.slopes <- fish.diversity %>%
  mutate(id = paste(longitude, latitude, sep = "_"),
         year = as.numeric(as.character(year))) %>%
  split(.$id) %>%
  purrr::map(~lm(fric_pa ~ year, data = .x)) %>% 
  purrr::map_df(broom::tidy, .id = 'id') %>%
  filter(term == 'year') %>%
  tidyr::separate(id, c("longitude", "latitude"), sep = "_") %>%
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude),
         lowerCI = estimate - 1.96 * std.error,
         upperCI = estimate + 1.96 * std.error)
slopeLims <- round(max(abs(range(r.slopes$estimate))) * c(-1, 1), 3)
plot9 <- plotFun(r.slopes, variable = "estimate", size = 0.5) +
  scale_color_gradientn(colours = rev(brewer.pal(11, "RdYlBu")),
                        limits = slopeLims, name = "Variation") ; plot9

fish.diversity_mean <- cbind(fish.diversity_mean,
                        "Trend fric_pa" = r.slopes$estimate
)

# saveRDS(fish.diversity_mean, file = "./data/Biodiversity metrics/fish.diversity_mean.rds")

data <- fish.diversity_mean

plotFun <- function(data, variable, years = NULL, size = 0.5){
  if(!is.null(years)){
    data <- data %>%
      filter(year %in% years)
  }
  p <- ggplot(data, aes(x = longitude, y = latitude, color = .data[[variable]])) +
    geom_point(size = size) + 
    borders(fill = "grey", colour = "black", size = 0.25) + 
    #coord_quickmap(xlim = range(data$longitude), ylim = range(data$latitude)) +
    coord_quickmap(xlim = c(9.041542, 23.179892), ylim = c(53.59160 , 59.85796)) +
    scale_color_gradientn(colours = rev(brewer.pal(11, "RdYlBu"))) +
    #annotate("text", x = min(data$longitude), y = max(data$latitude), label = "A", hjust = 0, vjust = 1, size = 5) +
    annotate("text", x = 9.041542, y = 59.85796, label = "A", hjust = 0, vjust = 1, size = 5) +
    theme(
      legend.position = c(0.95, 0.15),
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_blank()
    )
  
  if(!is.null(years)){
    p <- p + facet_wrap(~ year)
  }
  return(p)
}

# A
plot1 <- plotFun(fish.diversity_mean, variable = "Richness", years = c()) ; plot1

# C
plot2 <- plotFun(fish.diversity_mean, variable = "Shannon", years = c()) ; plot2

# E
plot3 <- plotFun(fish.diversity_mean, variable = "funHill0", years = c()) ; plot3

# G
plot4 <- plotFun(fish.diversity_mean, variable = "fric_pa", years = c()) ; plot4

# I
plot5 <- plotFun(fish.diversity_mean, variable = "fdiv", years = c()) ; plot5

# B
plot6 <- plotFun(fish.diversity_mean, variable = "Trend Richness", years = c()) ; plot6

# D
plot7 <- plotFun(fish.diversity_mean, variable = "Trend Shannon", years = c()) ; plot7

# F
plot8 <- plotFun(fish.diversity_mean, variable = "Trend funHill0", years = c()) ; plot8

# H
plot9 <- plotFun(fish.diversity_mean, variable = "Trend fric_pa", years = c()) ; plot9

# J
plot10 <- plotFun(fish.diversity_mean, variable = "Trend fdiv", years = c()) ; plot10
