################################################################################
##### Script 6 : Make predictions ----------------------------------------------
##### Tom Jimenez, 16-11-2026 --------------------------------------------------
################################################################################

remove(list=ls())

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")
grid.directory  <- file.path(data.directory, "grid_Copernicus_Baltic_1993-2021")
out.directory   <- file.path(data.directory, "EpredY")
dir.create(out.directory, showWarnings = FALSE, recursive = TRUE)

library(dplyr)
library(sp)
library(Hmsc)
set.seed(1)

nChains = 4
samples = 500
thin = 500
filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin),"_occurrence"))
#filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin),"_biomass"))
load(filename)

# keep 1 posterior sample out of 10 -> 200 samples
# 2000 samples in total = 4 chains x 500
# Larger value = faster and slightly noisier maps.
thin_post <- 10 # 1 if you want all samples
# I took 1 for the results in the paper

# Transformation applied to the biomass before fitting the biomass model: Y = log_base-log(B + offset)
#   log(B + 1)   -> log_base = exp(1), log_offset = 1   (this model)
#   log10(B + 1) -> log_base = 10,      log_offset = 1
# Absences (B = 0) were set to NA, so the model describes biomass given presence.
log_base   <- exp(1)
log_offset <- 1

# Upper bound on the predicted biomass of each species, as a multiple of the
# MAXIMUM biomass observed for that species in the fitting data.
# WHY: biomass = exp(linear predictor). Where the covariates of a grid cell are
# outside the range the species was observed in (polynomial terms!), the linear
# predictor can reach 100+ and exp() gives absurd values (1e50, Inf). The
# probit occurrence model cannot do this (probabilities are bounded by 1).
#   1   = never above the maximum observed biomass (recommended)
#   2   = more permissive;   Inf = no bound (NOT recommended)
cap_factor <- 1

years  <- 1993:2021
# Months present in the fitting data. There are no observations in April, May
# and June (no "Spring" level in the model), so these months cannot be predicted.
months <- c(1, 2, 3, 7, 8, 9, 10, 11, 12)

# Season of each month, indexed by month number (1-12).
# Must be the same month -> season mapping as the one used to fit the model.
season_of_month <- c("Winter", "Winter", "Winter",   # Jan, Feb, Mar
                     NA,       NA,       NA,         # Apr, May, Jun (not in model)
                     "Summer", "Summer", "Summer",   # Jul, Aug, Sep
                     "Autumn", "Autumn", "Autumn")   # Oct, Nov, Dec

# File where the spatial latent factors of the grid are cached 
# so that the slow step is not repeated if the script is re-run
# DELETE this file if you change the model, the grid or thin_post
eta_file <- file.path(data.directory,
                      paste0("eta_cell_meanfield_", floor(nChains * samples / thin_post), ".rds"))
#eta_file <- file.path(data.directory,
#                      paste0("eta_cell_biomass_meanfield_", floor(nChains * samples / thin_post), ".rds"))

# Fitting data
data <- readRDS("./data/fish_survey_Baltic_sea_v2.rds")

# Calculate the longitude and latitude limits in `data`.
min_lon_data <- min(data$lon, na.rm = TRUE) - 1
max_lon_data <- max(data$lon, na.rm = TRUE) + 1
min_lat_data <- min(data$lat, na.rm = TRUE) - 1
max_lat_data <- max(data$lat, na.rm = TRUE) + 1

# ---- Sanity checks on the model ----------------------------------------------
# The code below assumes that random level 1 = "cell" and level 2 = "Year"
stopifnot(all(m$rLNames == c("cell", "Year")))

# Seasons known by the model (here: "Autumn" "Summer" "Winter")
model_seasons <- levels(m$XData$season)
stopifnot(all(season_of_month[months] %in% model_seasons))

# Every year to predict must be one of the Year units estimated by the model
# (so the estimated year effect is re-used, not interpolated)
stopifnot(all(as.character(years) %in% levels(m$dfPi$Year)))

# The spatial level of the model is stored as an sp::SpatialPoints object
stopifnot(inherits(m$rL$cell$s, "SpatialPoints"))

# ---- Occurrence model
# The model must be a probit model for every species (family 2 = probit)
stopifnot(all(m$distr[, 1] == 2))

# ---- Biomass model
# The model must be a NORMAL model for every species (family 1 = normal)
stopifnot(all(m$distr[, 1] == 1))

# Scaling of Y used when fitting: row 1 = mean, row 2 = sd of each species
# (0 and 1 if Y was not scaled). Used to go back to the scale of the data.
y_mean <- m$YScalePar[1, ]
y_sd   <- m$YScalePar[2, ]

# Maximum observed biomass of each species (original units), from the fitting
# data. m$Y is the data as given to Hmsc (already on the log(B + 1) scale and
# NOT standardised, even if YScale = TRUE). Used as an upper bound in the
# predictions (see cap_factor).
B_cap <- cap_factor * (log_base^apply(m$Y, 2, max, na.rm = TRUE) - log_offset)
stopifnot(length(B_cap) == m$ns, all(is.finite(B_cap) | is.infinite(cap_factor)))
# ------------------------------------------------------------------------------

# Path of the grid file of a given year and month
grid_file <- function(year, month) {
  file.path(grid.directory,
            paste0("grid_Copernicus_Baltic_", year, "_", sprintf("%02d", month), ".rds"))
}

# Read a grid file
read_grid <- function(year, month) {
  tryCatch(readRDS(grid_file(year, month)),
           error = function(e) {
             message("Cannot read: ", grid_file(year, month))
             NULL
           })
}

# Design matrix of the covariates, built with the model's own formula.
# IMPORTANT: Beta in m$postList is on the ORIGINAL scale of the covariates,
# so the design matrix must NOT be standardised.
# Columns are: (Intercept), poly(T_floor)1, poly(T_floor)2, ..., seasonWinter.
# This is safe because the formula uses poly(..., raw = TRUE).
make_X <- function(XD) model.matrix(m$XFormula, XD)

# Check that make_X reproduces exactly the model's own design matrix
stopifnot(all(colnames(make_X(m$XData)) == m$covNames))
stopifnot(max(abs(unname(make_X(m$XData)) - unname(m$X))) == 0)

# The Copernicus grid has the same cells every month, so cell selection and the
# spatial effect only need to be computed once. Use the first readable file as
# the reference grid.
raw0 <- NULL
for (y in years) {
  for (mo in months) {
    raw0 <- read_grid(y, mo)
    if (!is.null(raw0)) break
  }
  if (!is.null(raw0)) break
}
stopifnot(!is.null(raw0))

# Keep cells inside the study-area bounding box
in_box <- raw0$lon >= min_lon_data & raw0$lon <= max_lon_data &
  raw0$lat >= min_lat_data & raw0$lat <= max_lat_data

# Discard cells too far from observations (TRUE = too far)
too_far <- mgcv::exclude.too.far(raw0$lon[in_box], raw0$lat[in_box],
                                 data$lon, data$lat, dist = 0.01)

# `sel` = row indices (in the raw grid file) of the cells that are predicted.
# The same indices are applied to every monthly file.
sel   <- which(in_box)[!too_far]
grid0 <- raw0[sel, ]
n0    <- nrow(grid0)           # number of cells to predict
cat("Number of grid cells to predict:", n0, "\n")

# Save the coordinates of the predicted cells so that the maps can be drawn
# later: row_id matches the row names of every EpredY_*.rds file.
saveRDS(data.frame(row_id = sel, lon = grid0$lon, lat = grid0$lat),
        file.path(out.directory, "grid_coordinates.rds"))

# Spatial random level
# We start from the model's own spatial level (m$rL$cell) so that it keeps the
# NNGP settings, the number of neighbours, the longlat option and, above all,
# the alpha grid (alphapw) that the posterior alpha indices refer to.
# Then we append the new grid cells to its coordinates and unit names.

cell_ids <- paste0("g", seq_len(n0)) # unique names for the new cells

xyOld <- m$rL$cell$s # 47 fitted cells (SpatialPoints)
xyNew <- SpatialPoints(
  matrix(c(grid0$lon, grid0$lat), ncol = 2,
         dimnames = list(cell_ids, colnames(coordinates(xyOld)))),
  proj4string = xyOld@proj4string) # same projection (longlat WGS84)

rL.cell    <- m$rL$cell
rL.cell$s  <- rbind(xyOld, xyNew) # fitted cells first, then new cells
rL.cell$pi <- factor(c(as.character(m$rL$cell$pi), cell_ids))

# Spatial latent factor (eta) of the whole grid, computed once
# Posterior samples used for the predictions (pooled over the 4 chains)
post <- poolMcmcChains(m$postList, thin = thin_post)
cat("Number of posterior samples used:", length(post), "\n")

# Predict the spatial latent factors of the new cells, conditional on the
# estimated factors of the 47 fitted cells.
#   - predictMeanField = TRUE: fast vectorised conditional prediction (no loop
#     over cells). It was checked against the NNGP prediction and the
#     difference was the same size as the NNGP Monte Carlo noise.
#   - The NNGP branch of Hmsc would take many hours for 39,000 cells.
get_eta <- function(post, idx = seq_len(n0), ...) {
  Hmsc:::predictLatentFactor(
    unitsPred = cell_ids[idx],                 # new cells to predict
    units     = levels(m$dfPi$cell),           # fitted cells
    postEta   = lapply(post, function(p) p$Eta[[1]]),     # level 1 = cell
    postAlpha = lapply(post, function(p) p$Alpha[[1]]),
    rL        = rL.cell,
    ...)
}

if (file.exists(eta_file)) {
  # Re-use the cached result of a previous run
  cat("Loading cached spatial factors:", eta_file, "\n")
  eta_cell <- readRDS(eta_file)
} else {
  cat("Computing spatial factors for the whole grid (a few minutes)...\n")
  print(system.time(eta_cell <- get_eta(post, predictMeanField = TRUE)))
  saveRDS(eta_cell, eta_file)
}

# eta_cell is a list (one element per posterior sample) of matrices
# (n0 cells x number of factors). Check that it matches the current setup.
stopifnot(length(eta_cell) == length(post),
          nrow(eta_cell[[1]]) == n0)

# ---- Prediction functions
# ---- Function for the occurrence model
# Posterior mean occurrence probability for a set of cells in one year.
#   X       : design matrix (cells x covariates), NOT standardised
#   etaCell : list (one per posterior sample) of spatial factors (cells x nf)
#   post    : list of posterior samples
#   iy      : position of the year in the model's Year units
# For each posterior sample s, the probit linear predictor is
#   L = X Beta + Eta_cell Lambda_cell + Eta_year Lambda_year
# and the occurrence probability is pnorm(L). The result is averaged over the
# samples.
pred_core <- function(X, etaCell, post, iy) {
  E <- matrix(0, nrow(X), m$ns)
  for (s in seq_along(post)) {
    p <- post[[s]]
    # Fixed effects + spatial random effect (level 1 = cell)
    L <- X %*% p$Beta + etaCell[[s]] %*% p$Lambda[[1]]
    # Year random effect (level 2): one value per species, the same for all
    # cells of the month, added to every row of L
    L <- L + rep(as.numeric(p$Eta[[2]][iy, , drop = FALSE] %*% p$Lambda[[2]]),
                 each = nrow(L))
    E <- E + pnorm(L)
  }
  E / length(post)
}

# ---- Function for the biomass model
# Posterior mean biomass (given presence) for a set of cells in ONE year.
#   X       : design matrix (cells x covariates), NOT standardised
#   etaCell : list (one per posterior sample) of spatial factors (cells x nf)
#   post    : list of posterior samples
#   iy      : position of the year in the model's Year units
# For each posterior sample s, the linear predictor of the normal model is
#   L = X Beta + Eta_cell Lambda_cell + Eta_year Lambda_year
# which is the expected log(biomass + 1) (on the scale of Y used in the fit).
# Steps for each sample:
#   1. back to the scale of the transformed data: L * sd + mean
#      (does nothing if Y was not scaled: sd = 1, mean = 0)
#   2. residual variance of the data scale: sigma * sd^2
#      (post[[s]]$sigma is a VARIANCE, one value per species)
#   3. mean of a log-normal variable, in the ORIGINAL units of the data:
#        B = exp(k * L + k^2 * sigma / 2) - log_offset,   k = log(log_base)
#      (the "+ sigma/2" term is what makes B the MEAN and not the median;
#       the offset of the log(B + 1) transformation is subtracted at the end).
#      A predicted B slightly below 0 is not possible for a biomass: it is set
#      to 0 (this only happens for tiny biomasses, where L is close to 0).
#   4. B is capped at cap_factor x the maximum observed biomass of the species
#      (protects against exp() blowing up when the covariates are outside the
#      range of the data).
# The result is averaged over the samples.
pred_core <- function(X, etaCell, post, iy) {
  k <- log(log_base)
  B <- matrix(0, nrow(X), m$ns)
  for (s in seq_along(post)) {
    p <- post[[s]]
    # Fixed effects + spatial random effect (level 1 = cell)
    L <- X %*% p$Beta + etaCell[[s]] %*% p$Lambda[[1]]
    # Year random effect (level 2): one value per species, the same for all
    # cells of the month, added to every row of L
    L <- L + rep(as.numeric(p$Eta[[2]][iy, , drop = FALSE] %*% p$Lambda[[2]]),
                 each = nrow(L))
    # 1. back to the scale of the transformed data
    L <- L * rep(y_sd, each = nrow(L)) + rep(y_mean, each = nrow(L))
    # 2. residual variance on the same scale (one value per species)
    sig <- p$sigma * y_sd^2
    # 3. mean biomass given presence, original units (offset removed, >= 0)
    #    (bounded above by the maximum observed biomass of each species)
    B <- B + pmin(pmax(exp(k * L + rep(k^2 * sig / 2, each = nrow(L))) -
                         log_offset, 0),
                  rep(B_cap, each = nrow(L)))
  }
  B / length(post)
}

# Predict one month. Returns a (cells x species) matrix, or NULL if the grid
# file is unreadable or differs from the reference grid.
predict_month <- function(year, month) {
  
  g <- read_grid(year, month)
  if (is.null(g)) return(NULL)
  
  # The shared spatial effect is only valid if the grid cells are identical
  if (!identical(g[, c("lon", "lat")], raw0[, c("lon", "lat")])) {
    message("Grid differs from the reference grid for ", year, "-", month,
            ": skipped")
    return(NULL)
  }
  g <- g[sel, ] # same cells as the reference
  
  # Covariates of the month. Column names must match the model's formula.
  # `season` must have exactly the levels of the model.
  XD <- data.frame(
    Sal_floor = g$Sal_floor, T_floor = g$T_floor, O2_floor = g$O2_floor,
    depth     = g$depth,     chl     = g$chl,
    season    = factor(rep(season_of_month[month], nrow(g)),
                       levels = model_seasons))
  
  # model.matrix drops rows with missing covariates, so remove them explicitly
  # to keep rows aligned with the spatial factors
  keep <- complete.cases(XD)
  if (!any(keep)) return(NULL)
  
  E <- pred_core(make_X(XD[keep, ]),
                 lapply(eta_cell, function(e) e[keep, , drop = FALSE]),
                 post,
                 match(as.character(year), levels(m$dfPi$Year)))
  
  # Row names = row index in the raw grid file; column names = species
  dimnames(E) <- list(sel[keep], m$spNames)
  E
}

# Loop over years and months
for (year in years) {
  for (month in months) {
    
    out_file <- file.path(out.directory,
                          paste0("EpredY_", year, "_", sprintf("%02d", month), ".rds"))
    #out_file <- file.path(out.directory,
    #                      paste0("EpredY_biomass_", year, "_", sprintf("%02d", month), ".rds"))
    if (file.exists(out_file)) next
    
    EpredY <- predict_month(year, month)
    if (is.null(EpredY)) next
    
    saveRDS(EpredY, out_file)
    cat("Predictions for year", year, "month", sprintf("%02d", month), "finished.\n")
    remove(EpredY)
    gc()
  }
}

# Check which months are missing
todo <- expand.grid(year = years, month = months)
todo$file <- file.path(out.directory,
                       paste0("EpredY_", todo$year, "_", sprintf("%02d", todo$month), ".rds"))
#todo$file <- file.path(out.directory,
#                       paste0("EpredY_biomass_", todo$year, "_", sprintf("%02d", todo$month), ".rds"))
missing_months <- todo[!file.exists(todo$file), c("year", "month")]

cat("\nMonths computed:", sum(file.exists(todo$file)), "out of", nrow(todo), "\n")
if (nrow(missing_months) > 0) {
  cat("Missing months (unreadable/truncated grid files, to be copied again):\n")
  print(missing_months)
}

# ---- Combine all EpredY ------------------------------------------------------

years  <- 1993:2021
months <- c(1, 2, 3, 7, 8, 9, 10, 11, 12)

# Initialise an empty DataFrame to store all rows
combined_data <- data.frame()

# Loop to load each EpredY file and add it to the DataFrame
for (year in years) {
  for (month in months) {
    
    # Format month with leading zero (01, 02, ..., 12)
    month_str <- sprintf("%02d", month)
    
    # Load EpredY
    EpredY <- readRDS(paste0("./data/EpredY/EpredY_", year, "_", month_str, ".rds"))
    #EpredY <- readRDS(paste0("./data/EpredY/EpredY_biomass_", year, "_", month_str, ".rds"))
    EpredY <- as.data.frame(EpredY)
    
    # Add a ‘year’ and a ‘month’ columns to identify each observation
    EpredY$Year <- year
    EpredY$Month <- month
    
    # Add lines from the current file to the combined DataFrame
    combined_data <- rbind(combined_data, EpredY)
    cat("Data for the year", year, "month", month_str, "added.\n")
  }
}

grid_coordinates <- readRDS(file.path(out.directory, "grid_coordinates.rds"))

expanded_grid <- grid_coordinates[rep(1:nrow(grid_coordinates), times = 261), ] # 261, as 29 years multiplied by 9 months
expanded_grid <- expanded_grid[1:nrow(combined_data), ]

combined_data$lon <- expanded_grid$lon
combined_data$lat <- expanded_grid$lat

combined_data <- combined_data %>%
  select(lon, lat, Year, Month, everything())

names(combined_data)[names(combined_data) == "lon"] <- "Longitude"
names(combined_data)[names(combined_data) == "lat"] <- "Latitude"

saveRDS(combined_data, file = paste0("./data/EpredY/Combined_EpredY_1993_01-2021_12.rds"))
#saveRDS(combined_data, file = paste0("./data/EpredY/Combined_EpredY_biomass_1993_01-2021_12.rds"))
cat("The combined file has been saved.\n")

# Epredy Average per year
for (year in years) {
  
  month_str <- sprintf("%02d", months)
  paths <- paste0("./data/EpredY/EpredY_", year, "_", month_str, ".rds")
  #paths <- paste0("./data/EpredY/EpredY_biomass_", year, "_", month_str, ".rds")
  paths_exist <- paths[file.exists(paths)]
  
  if (length(paths_exist) == 0) {
    cat("No files found for year", year, "— skipping.\n")
    next
  }
  cat("Year", year, "— averaging", length(paths_exist), "months... ")
  
  EpredY_sum <- NULL
  
  for (p in paths_exist) {
    EpredY_month <- readRDS(p)
    if (is.null(EpredY_sum)) {
      EpredY_sum <- EpredY_month
    } else {
      EpredY_sum <- EpredY_sum + EpredY_month
    }
    remove(EpredY_month)
  }
  
  # Mean
  EpredY_annual <- EpredY_sum / length(paths_exist)
  
  # Save
  saveRDS(EpredY_annual, paste0("./data/EpredY/EpredY_", year, ".rds"))
  #saveRDS(EpredY_annual, paste0("./data/EpredY/EpredY_biomass_", year, ".rds"))
  
  cat("EpredY for the year", year, "saved.\n")
  remove(EpredY_sum, EpredY_annual)
  gc()
}

cat("\nDone —", length(years), "annual EpredY computed.\n")

# Combine
combined_data <- data.frame()

for (year in years) {
  
  # Load EpredY
  EpredY <- readRDS(paste0("./data/EpredY/EpredY_", year,".rds"))
  #EpredY <- readRDS(paste0("./data/EpredY/EpredY_biomass_", year,".rds"))
  EpredY <- as.data.frame(EpredY)
  
  # Add a ‘year’ column to identify each observation
  EpredY$Year <- year
  
  # Add lines from the current file to the combined DataFrame
  combined_data <- rbind(combined_data, EpredY)
  cat("Data for the year", year,"added.\n")
}

grid_coordinates <- readRDS(file.path(out.directory, "grid_coordinates.rds"))

expanded_grid <- grid_coordinates[rep(1:nrow(grid_coordinates), times = 29), ] # 29 years
expanded_grid <- expanded_grid[1:nrow(combined_data), ]

combined_data$lon <- expanded_grid$lon
combined_data$lat <- expanded_grid$lat

combined_data <- combined_data %>%
  select(lon, lat, Year, everything())

names(combined_data)[names(combined_data) == "lon"] <- "Longitude"
names(combined_data)[names(combined_data) == "lat"] <- "Latitude"

saveRDS(combined_data, file = paste0("./data/EpredY/Combined_EpredY_1993-2021.rds"))
#saveRDS(combined_data, file = paste0("./data/EpredY/Combined_EpredY_biomass_1993-2021.rds"))
cat("The combined file has been saved.\n")

# Epredy Average
# Initialize a running sum DataFrame and a count DataFrame
sum_EpredY <- NULL
count_EpredY <- NULL

# Loop to load each EpredY file and update the sum and count
for (year in years) {
  # Load the EpredY file for the current year
  EpredY <- readRDS(paste0("./data/EpredY/EpredY_", year, ".rds"))
  #EpredY <- readRDS(paste0("./data/EpredY/EpredY_biomass_", year, ".rds"))
  EpredY <- as.as.matrix(EpredY)
  
  ok <- !is.na(EpredY)          # valid values of this year
  EpredY[!ok] <- 0              # NA contributes 0 to the sum
  
  # If it's the first year, initialize the sum and count DataFrames
  if (is.null(sum_EpredY)) {
    sum_EpredY <- EpredY
    count_EpredY <- ok + 0L     # logical matrix -> integer matrix
  } else {
    # Add the current year to the sum, ignoring NA values
    # Stop if the rows/columns are not the same cells/species as the first year
    stopifnot(identical(dim(EpredY), dim(sum_EpredY)),
              identical(rownames(EpredY), rownames(sum_EpredY)),
              identical(colnames(EpredY), colnames(sum_EpredY)))
    sum_EpredY   <- sum_EpredY + EpredY
    count_EpredY <- count_EpredY + ok
  }
  cat("Data for year", year, "processed.\n")
}

# Mean over the valid years; NA only where no year has a valid value
EpredY_mean <- sum_EpredY / count_EpredY
EpredY_mean[count_EpredY == 0] <- NA

# Convert to DataFrame for easier handling
EpredY_mean <- as.data.frame(EpredY_mean)

EpredY_mean$lon <- grid_coordinates$lon
EpredY_mean$lat <- grid_coordinates$lat

EpredY_mean <- EpredY_mean %>%
  select(lon, lat, everything())

names(EpredY_mean)[names(EpredY_mean) == "lon"] <- "Longitude"
names(EpredY_mean)[names(EpredY_mean) == "lat"] <- "Latitude"

saveRDS(EpredY_mean, file = paste0("./data/EpredY/EpredY_mean_1993-2021.rds"))
#saveRDS(EpredY_mean, file = paste0("./data/EpredY/EpredY_biomass_mean_1993-2021.rds"))
cat("The final file has been saved.\n")
