################################################################################
##### Script 7 : Biodiversity metrics and trends -------------------------------
##### Tom Jimenez, 16-11-2026 --------------------------------------------------
################################################################################

remove(list=ls())

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

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
biomass <- readRDS("./data/EpredY/Combined_EpredY_biomass_1993-2021.rds")
#p.occurrence <- readRDS("./data/EpredY/EpredY_mean_1993-2021.rds")
#biomass <- readRDS("./data/EpredY/EpredY_biomass_mean_1993-2021.rds")

p.occurrence[1:5, 1:10] %>% kable()
biomass[1:5, 1:10] %>% kable()

# Are all columns in the same order?  

table(names(p.occurrence) == names(biomass)) #TRUE

all(biomass %>% select(Year, Longitude, Latitude) == p.occurrence %>% select(Year, Longitude, Latitude)) # TRUE

# Create vector with species names

species <- names(p.occurrence)[!names(p.occurrence) %in% c("Year", "Longitude", "Latitude")]
species[1:10]

##### Taxonomic diversity ######################################################
##### Taxonomic richness -------------------------------------------------------

# Expected species richness = sum of the occurrence probabilities of all species in each row
# Keep the time and location variables, add the richness column
richness <- cbind(p.occurrence %>% select(-all_of(species)),
                  Richness = rowSums(p.occurrence %>% select(all_of(species))))

range(richness$Richness)

##### Shannon & evenness -------------------------------------------------------

# Unconditional expected biomass = P(presence) x biomass given presence
biomass[, species] <- biomass[, species] * p.occurrence[, species]

shannon <- diversity(biomass[, species], index = "shannon")
range(shannon)
shannon[1:10]

# Compute evenness
# Pielou over the full species pool (= Shannon rescaled to 0-1)

evenness <- shannon / log(length(species))

# Effective number of species (same scale as richness)
effective_sp <- exp(shannon)
range(effective_sp)

# Merge with previous biodiversity metrics and visualize

all(biomass %>% select(Year, Longitude, Latitude) == richness %>% select(Year, Longitude, Latitude)) # TRUE

fish.diversity <- cbind(richness,
                        EffectiveSpecies = effective_sp,
                        Shannon = shannon,
                        Evenness = evenness) %>%
  mutate(Year = as.numeric(as.character(Year)))

ggplot(fish.diversity, aes(x = Richness, y = EffectiveSpecies)) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4))

ggplot(fish.diversity, aes(x = Richness, y = EffectiveSpecies)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), color = "red")

##### Functional diversity #####################################################

# The computation of functional diversity is based on the mFD package
# Essentially there are two ways to obtain different functional biodiversity metrics:
# 1 - Using the Hill numbers
# 2 - Using a multidimensional space based on a PCoA of species traits

##### Hill numbers -------------------------------------------------------------

# Load species traits

load("./data/TrData.RData")
TrData[1:5, 1:8] %>% kable()

# Explore possible correlations among trais, here we set a threshold of 0.7

TrData %>%
  select_if(is.numeric) %>%
  cor() %>%
  as.data.frame() %>%
  mutate_all(.funs = function(x) round(x, 2)) %>%
  mutate_all(.funs = function(x) ifelse(abs(x) < 0.7, NA, x)) %>%
  kable()

# Make sure that species names (i.e., taxon) is set as rownames

TrData <- as.data.frame(TrData)
rownames(TrData) <- TrData$taxon
TrData$taxon <- NULL
TrData[1:5, 1:7] %>% kable()

# Retrieve trait names and type (i.e., ordinal, categorical continuous, fuzzy) into a table

tr_cat <- tibble(trait_name  = names(TrData),
                 trait_type = c(rep("N", times = 3), # Nominal traits (factor variable)
                                rep("Q", times = 4))) # Quantitative traits (numeric values)
head(tr_cat) %>% kable()

# Compute species distances

sp_dist <- funct.dist(sp_tr = TrData,
                      tr_cat = tr_cat, 
                      # Since there are categorical traits, we compute gower distance
                      metric = "gower", 
                      # If only continuous traits, one can use euclidean distance instead
                      # metric = "euclidean",
                      scale_euclid = "scale_center")

as.matrix(sp_dist)[1:5, 1:10] %>% kable()

# Now one can use these functional distances to compute a series of functional diversity metrics (Hill numbers)

rownames(biomass) <- NULL
biomass.assemblage <- biomass[, species] %>%
  mutate(ID = paste("community", 1:nrow(biomass), sep = "_")) %>%
  column_to_rownames("ID") %>%
  as.matrix()

biomass.assemblage[1:5, 1:10] %>% kable()

hill_test <- alpha.fd.hill(asb_sp_w = biomass.assemblage[1:1000,],
                           sp_dist = sp_dist,
                           tau = "mean",
                           q = c(1, 2))
head(hill_test$asb_FD_Hill) %>% kable()

hill12 <- alpha.fd.hill(asb_sp_w = biomass.assemblage,
                        sp_dist = sp_dist,
                        tau = "mean",
                        q = c(1, 2)) # Shannon-like & # Simpson-like

# Make sure that the data used when computing Hill numbers has same order than other biodiversity metrics, and merge

fish.diversity <- cbind(fish.diversity,
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

# mad with 4 dimensions has the smallest value => this are the best settings to represent the community from a functional perspective
# We can also visualize the quality of the selected functional spaces

quality.fspaces.plot(fspaces_quality = quality_fspaces,
                     quality_metric = "mad",
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

# As described previously, the PCoA with mad and 4 dimensions is the setup that best describes the community from a functional perspective. Let’s proceed with that

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

tr_faxes_Q <- traits.faxes.cor(sp_tr = TrData[, tr_cat$trait_name[tr_cat$trait_type == "Q"]], # Only QUANTITATIVE traits
                               sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3")], # Display only relationships with the first 3 axis (PCs)
                               plot = TRUE) # return plot
tr_faxes_Q

tr_faxes_N <- traits.faxes.cor(sp_tr = TrData[, tr_cat$trait_name[tr_cat$trait_type == "N"]], # Only CATECORICAL traits
                               sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3")],
                               plot = TRUE)
tr_faxes_N

# Finally we can compute the functional diversity metrics based on the PCoA. A reminder that even though the PCoA has 8 PCs, we only use 5 to slow down computation

biomass.assemblage[1:5, 1:10] %>% kable()

# Test with few observations first

test_fundiv <- alpha.fd.multidim(
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5")],
  asb_sp_w = biomass.assemblage[1:100, ],
  # Define which functional diversity metrics are to be computed
  ind_vect = c("feve", "fdis", "fdiv"),
  scaling = TRUE,
  check_input = TRUE,
  details_returned = TRUE)

head(test_fundiv$functional_diversity_indices) %>% kable()

# Running the whole set of communities (i.e., rows) 

alpha_fd_indices <- alpha.fd.multidim(
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5")],
  asb_sp_w = biomass.assemblage,
  # Define which functional diversity metrics are to be computed
  # Note that more indices will slow down the process
  ind_vect = c("feve", "fdis", "fdiv"),
  scaling = TRUE,
  check_input = TRUE,
  details_returned = TRUE)

head(alpha_fd_indices$functional_diversity_indices) %>% kable()

# We can finally merge the multivariate functional diversity metrics with all other

fish.diversity <- cbind(fish.diversity,
                        alpha_fd_indices$functional_diversity_indices)

fish.diversity[1:5, 1:10] %>% kable()

# Another way of calculating indices by subset
# Calculation of hill12 by subsets 

n_rows <- nrow(biomass.assemblage)
chunk_size <- 100000      
output_dir <- "./data/Biodiversity_metrics/output_hill12_rds" 
dir.create(output_dir, showWarnings = FALSE) 

for (start_row in seq(1, n_rows, by = chunk_size)) {
  end_row <- min(start_row + chunk_size - 1, n_rows)
  message(paste("Processing of lines", start_row, "to", end_row))
  hill_result <- alpha.fd.hill(
    asb_sp_w = biomass.assemblage[start_row:end_row, ],
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

saveRDS(hill_combined, "./data/Biodiversity_metrics/Hill12_combined.rds")

hill12 <- read_rds("./data/Biodiversity_metrics/Hill12_combined.rds")
fish.diversity <- cbind(fish.diversity,
                        funHill1 = hill12$FD_q1,
                        funHill2 = hill12$FD_q2)

# Calculation of alpha_fd_indices by subsets 

n_rows <- nrow(biomass.assemblage)
chunk_size <- 10000      
output_dir <- "./data/Biodiversity_metrics/output_alpha_fd_indices_rds" 
dir.create(output_dir, showWarnings = FALSE) 

for (start_row in seq(1, n_rows, by = chunk_size)) {
  end_row <- min(start_row + chunk_size - 1, n_rows)
  message(paste("Processing of lines", start_row, "to", end_row))
  fundiv_result <- alpha.fd.multidim(
    sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5")],
    asb_sp_w = biomass.assemblage[start_row:end_row, ],
    ind_vect = c("feve", "fdis", "fdiv"),
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

saveRDS(fundiv_combined, "./data/Biodiversity_metrics/alpha_fd_indices_combined.rds")

alpha_fd_indices <- read_rds("./data/Biodiversity_metrics/alpha_fd_indices_combined.rds")
fish.diversity <- cbind(fish.diversity,
                        feve = alpha_fd_indices$feve,
                        fdis = alpha_fd_indices$fdis,
                        fdiv = alpha_fd_indices$fdiv,
                        fide_PC1 = alpha_fd_indices$fide_PC1,
                        fide_PC2 = alpha_fd_indices$fide_PC2,
                        fide_PC3 = alpha_fd_indices$fide_PC3,
                        fide_PC4 = alpha_fd_indices$fide_PC4,
                        fide_PC5 = alpha_fd_indices$fide_PC5)

#fish.diversity$community <- NULL

# saveRDS(fish.diversity, file = "./data/Biodiversity_metrics/fish.diversity.rds")
# fish.diversity <- read_rds("./data/Biodiversity_metrics/fish.diversity.rds")
# saveRDS(biomass.assemblage, file = "./data/Biodiversity_metrics/biomass.assemblage.rds")
# biomass.assemblage <- read_rds("./data/Biodiversity_metrics/biomass.assemblage.rds")

# Plots

fish.diversity <- read_rds("./data/Biodiversity_metrics/fish.diversity.rds")
fish.diversity_mean <- read_rds("./data/Biodiversity_metrics/fish.diversity_mean.rds")

plotFun <- function(data, variable, years = NULL, size = 0.5){
  if(!is.null(years)){
    data <- data %>%
      filter(Year %in% years)
  }
  p <- ggplot(data, aes(x = Longitude, y = Latitude, color = .data[[variable]])) +
    geom_point(size = size) + 
    borders(fill = "grey", colour = "black", size = 0.25) + 
    coord_quickmap(xlim = range(data$Longitude), ylim = range(data$Latitude)) +
    scale_color_gradientn(colours = rev(brewer.pal(11, "RdYlBu"))) +
    annotate("text", x = min(data$Longitude), y = max(data$Latitude), label = "A", hjust = 0, vjust = 1, size = 5) +
    theme(
      legend.position = c(0.95, 0.15),
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_blank()  
    )
  
  if(!is.null(years)){
    p <- p + facet_wrap(~ Year)
  }
  return(p)
}

# Compute linear trend over time
r.slopes <- fish.diversity %>%
  mutate(id = paste(Longitude, Latitude, sep = "_"),
         Year = as.numeric(as.character(Year))) %>%
  split(.$id) %>%
  purrr::map(~lm(Richness ~ Year, data = .x)) %>% 
  purrr::map_df(broom::tidy, .id = 'id') %>%
  filter(term == 'Year') %>%
  tidyr::separate(id, c("Longitude", "Latitude"), sep = "_") %>%
  mutate(Longitude = as.numeric(Longitude),
         Latitude = as.numeric(Latitude),
         lowerCI = estimate - 1.96 * std.error,
         upperCI = estimate + 1.96 * std.error)
slopeLims <- round(max(abs(range(r.slopes$estimate))) * c(-1, 1), 3)
plot <- plotFun(r.slopes, variable = "estimate", size = 0.5) +
  scale_color_gradientn(colours = rev(brewer.pal(11, "RdYlBu")),
                        limits = slopeLims, name = "Variation") ; plot

fish.diversity_mean <- cbind(fish.diversity_mean,
                        "Trend richness" = r.slopes$estimate
)

#saveRDS(fish.diversity_mean, file = "./data/Biodiversity_metrics/fish.diversity_mean.rds")

plotFun <- function(data, variable, label = NULL, years = NULL, size = 0.5){
  if(!is.null(years)){
    data <- data %>%
      filter(Year %in% years)
  }
  p <- ggplot(data, aes(x = Longitude, y = Latitude, color = .data[[variable]])) +
    geom_point(size = size) + 
    borders(fill = "grey", colour = "black", size = 0.25) + 
    #coord_quickmap(xlim = range(data$Longitude), ylim = range(data$Latitude)) +
    coord_quickmap(xlim = c(9.041542, 23.179892), ylim = c(53.59160 , 59.85796)) +
    scale_color_gradientn(colours = rev(brewer.pal(11, "RdYlBu"))) +
    #annotate("text", x = min(data$Longitude), y = max(data$Latitude), label = "A", hjust = 0, vjust = 1, size = 5) +
    annotate("text", x = 9.041542, y = 59.85796, label = label, hjust = 0, vjust = 1, size = 5) +
    theme(
      legend.position = c(0.95, 0.15),
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_blank()
    )
  
  if(!is.null(years)){
    p <- p + facet_wrap(~ Year)
  }
  return(p)
}

# A
plotFun(fish.diversity_mean, variable = "Richness", label = "A", years = c())
# B
plotFun(fish.diversity_mean, variable = "Trend richness", label = "B", years = c())

# C
plotFun(fish.diversity_mean, variable = "Shannon", label = "C", years = c())
# D
plotFun(fish.diversity_mean, variable = "Trend Shannon", label = "D", years = c())

# E
plotFun(fish.diversity_mean, variable = "funHill1", label = "E", years = c())
# F
plotFun(fish.diversity_mean, variable = "Trend funHill1", label = "F", years = c())

# G
plotFun(fish.diversity_mean, variable = "feve", label = "G", years = c())
# H
plotFun(fish.diversity_mean, variable = "Trend feve", label = "H", years = c())

# I
plotFun(fish.diversity_mean, variable = "fdis", label = "I", years = c())
# J
plotFun(fish.diversity_mean, variable = "Trend fdis", label = "J", years = c())

# K
plotFun(fish.diversity_mean, variable = "fdiv", label = "K", years = c())
# L
plotFun(fish.diversity_mean, variable = "Trend fdiv", label = "L", years = c())

###############################################################################
# Combined figure: one row per biodiversity metric,
# columns = [metric name | Mean map | Trend over time map]
###############################################################################

library(ggplot2)
library(patchwork)      # install.packages("patchwork") if needed
library(dplyr)
library(grid)           # textGrob, gpar
library(RColorBrewer)   # brewer.pal (used in plotFun)

# ---- 1. Map function (unchanged) --------------------------------------------
# plotFun: legend inside the panel + optional symmetric scale for trends
plotFun <- function(data, variable, label = NULL, years = NULL, size = 0.5,
                    symmetric = FALSE){
  if(!is.null(years)) data <- data %>% filter(Year %in% years)
  
  # Symmetric limits around 0 (so that yellow = no change)
  lim <- if (symmetric) {
    max(abs(range(data[[variable]], na.rm = TRUE))) * c(-1, 1)
  } else NULL
  
  p <- ggplot(data, aes(x = Longitude, y = Latitude, color = .data[[variable]])) +
    geom_point(size = size) +
    borders(fill = "grey", colour = "black", size = 0.25) +
    coord_quickmap(xlim = c(9.041542, 23.179892), ylim = c(53.59160, 59.85796)) +
    scale_color_gradientn(colours = rev(brewer.pal(11, "RdYlBu")), limits = lim) +
    annotate("text", x = 9.041542, y = 59.85796, label = label,
             hjust = 0, vjust = 1, size = 5) +
    theme(
      legend.position = c(0.99, 0.01),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_blank()
    )
  if(!is.null(years)) p <- p + facet_wrap(~ Year)
  p
}

# ---- 2. Table: one row per biodiversity metric -------------------------------
metrics <- tibble::tribble(
  ~row_label,                 ~mean_var,   ~trend_var,       ~lab_mean, ~lab_trend,
  "Taxonomic\nrichness",      "Richness",  "Trend richness", "A",       "B",
  "Shannon",                  "Shannon",   "Trend Shannon",  "C",       "D",
  "Functional\nHill (q = 1)", "funHill1",  "Trend funHill1", "E",       "F",
#  "Functional\nevenness",     "feve",      "Trend feve",     "G",       "H",
  "Functional\ndispersion",   "fdis",      "Trend fdis",     "G",       "H"
#  "Functional\ndivergence",   "fdiv",      "Trend fdiv",     "K",       "L"
)

# All variables must exist in the data
stopifnot(all(c(metrics$mean_var, metrics$trend_var) %in% names(fish.diversity_mean)))

# ---- 3. Helpers for the text panels ------------------------------------------
# Column header: horizontal text in its own panel
text_panel <- function(txt, size = 5, face = "bold") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = txt, size = size, fontface = face) +
    theme_void()
}

# Row label: rotated text drawn as a grob (not clipped by the panel)
label_panel <- function(txt, size = 11) {
  wrap_elements(full = textGrob(txt, rot = 90,
                                gp = gpar(fontface = "bold", fontsize = size,
                                          lineheight = 0.9)))
}

# ---- 4. Build the panels, row by row -----------------------------------------
# Layout: 3 columns = [row label | Mean | Trend over time]
panels <- list(plot_spacer(),                   # empty top-left corner
               text_panel("Mean"),              # column header 1
               text_panel("Trend over time"))   # column header 2

for (i in seq_len(nrow(metrics))) {
  panels <- c(panels, list(
    label_panel(metrics$row_label[i]),                         # metric name
    plotFun(fish.diversity_mean, variable = metrics$mean_var[i],
            label = metrics$lab_mean[i]),                      # Mean map
    plotFun(fish.diversity_mean, variable = metrics$trend_var[i],
            label = metrics$lab_trend[i], symmetric = TRUE)))  # Trend map
}

# ---- 5. Assemble ---------------------------------------------------------------
# widths : label column + two equal map columns
# heights: thin header row + one row per metric
final_plot <- wrap_plots(panels, ncol = 3,
                         widths  = c(0.1, 1, 1),
                         heights = c(0.06, rep(1, nrow(metrics))))

final_plot

# ---- 6. Save (adjust height if the maps look too stretched) -------------------
ggsave("./plots/fig_biodiversity_mean_trend.png", final_plot,
       width = 8, height = 12, dpi = 300)

dev.off()

# Correlations between metrics  ----

fish.diversity_mean %>%
  select(Richness, Shannon, funHill1, feve, fdis, fdiv) %>%
  cor(method = "spearman", use = "complete.obs") %>%
  round(2)

library(FactoMineR)
library(factoextra)

# ---- 1. PCA of the Mean maps (p1) ---------------------------------------------
selected_data <- fish.diversity_mean[, c("Richness", "Shannon", "funHill1", "fdis")]
colnames(selected_data) <- c("Taxonomic\nrichness", "Shannon",
                             "Functional\nHill (q = 1)", "Functional\ndispersion")

res_pca <- PCA(selected_data, graph = FALSE)

p1 <- fviz_pca_var(res_pca,
                   col.var = "contrib",
                   gradient.cols = c("blue", "red", "orange"),
                   repel = TRUE,
                   axes = c(1, 2)) +
  labs(title = NULL) +                        # remove the default title
  theme(legend.position = c(0.2, 0.4),
        legend.background = element_rect(fill = "white", color = "white")) +
  annotate("text", x = -Inf, y = Inf, label = "I",
           hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

# ---- 2. PCA of the Trend maps (p2) --------------------------------------------
selected_data <- fish.diversity_mean[, c("Trend richness", "Trend Shannon",
                                         "Trend funHill1", "Trend fdis")]
colnames(selected_data) <- c("Taxonomic\nrichness", "Shannon",
                             "Functional\nHill (q = 1)", "Functional\ndispersion")

res_pca <- PCA(selected_data, graph = FALSE)

p2 <- fviz_pca_var(res_pca,
                   col.var = "contrib",
                   gradient.cols = c("blue", "red", "orange"),
                   repel = TRUE,
                   axes = c(1, 2)) +
  labs(title = NULL) +
  theme(legend.position = c(0.2, 0.4),
        legend.background = element_rect(fill = "white", color = "white")) +
  annotate("text", x = -Inf, y = Inf, label = "J",
           hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

# ---- 3. Add the PCA row to the figure -----------------------------------------
# `panels` is the list built in the previous code (header + one row per metric).
# A new object is created so that re-running this block never adds the row twice.
panels_pca <- c(panels, list(label_panel("PCA"), p1, p2))

final_plot <- wrap_plots(panels_pca, ncol = 3,
                         widths  = c(0.1, 1, 1),
                         heights = c(0.06, rep(1, nrow(metrics)), 0.9))

final_plot

# ---- 4. Save: 4 map rows + 1 PCA row ------------------------------------------
ggsave("./plots/fig_biodiversity_mean_trend_pca.png", final_plot,
       width = 8, height = 16, dpi = 300)

dev.off()
