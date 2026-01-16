# Baltic Sea fish community HMSC model

This repository contains the scripts used in the manuscript Drivers and assembly processes shaping fish community composition and diversity: the Baltic Sea as a case study by Tom Jimenez, Marcel Montanyès, Federico Maioli, Baptiste Degueurce and Martin Lindegren.

## Scripts

- Script 1 - Descriptive figures and correlation matrix between parameters.
- Script 2 - Prepare all the data inputs for fitting the HMSC model, i.e., environmental data (X), community data (Y), trait data (T) and taxonomy (C). Also, define the model structure and define random effects. Finally, fit the presence/absence HMSC model.
- Script 2 bis - Same as script 2 but for the biomass HMSC model.
- Script 3 - Convergence evaluation.
- Script 4 - Compute model's explanatory and predictive power (4-fold crossvalidation), and exploration of model fit.
- Script 5 - Exploration of parameter estimates.
- Script 6 - Compute predictions and biodiversity metrics.

## Data sources

- The North Sea International Bottom Trawl Survey (NS-IBTS) and the Baltic International Trawl Survey (BITS) data are available from the DATRAS.
- Trait data were collected from available Beukhof et al. 2019 trait data base, supplemented with information from recent literature Coulon et al., 2023.
- Environmental data was retrieved from the physics and biochemistry reanalysis products of the Baltic Monitoring Forecasting Centre of Copernicus Marine Environment Monitoring Service (BALMFC CMEMS). The reanalysis products are based on the regionally downscaled Nemo-Nordic general circulation model (Hordoir et al. 2019) coupled with the SCOBI (Swedish Coastal and Ocean Biogeochemical) biogeochemistry model (Eilola et al. 2009; Almroth-Rosell et al. 2015).
