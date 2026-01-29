################################################################################
##### Script 3 : Evaluate convergence ------------------------------------------
################################################################################

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

library(Hmsc)
set.seed(1)
library(vioplot)
library(colorspace)

list.files(model.directory)

nChains = 4
samples = 250
thin = 100 
filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin),"_occurrence"))
#filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin),"_biomass"))
load(filename)

mpost = convertToCodaObject(m, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
tmp = mpost$Omega[[1]]
z = ncol(tmp[[1]])
sel = sample(z, size=100)
# Here we take the subset of species pairs. We loop over the 4 MCMC chains.
for(i in 1:length(tmp)){ 
  tmp[[i]] = tmp[[i]][,sel]
}
psrf.omega = gelman.diag(tmp,multivariate=FALSE)$psrf
psrf.gamma = gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf

ess.beta = effectiveSize(mpost$Beta)
ess.gamma = effectiveSize(mpost$Gamma)
ess.omega = effectiveSize(tmp)

mean(psrf.beta)
mean(ess.beta)
mean(psrf.gamma)
mean(ess.gamma)
mean(psrf.omega)
mean(ess.omega)

par(mfrow=c(1,3))
hist(psrf.beta, xlab = "psrf (beta)")
hist(psrf.gamma, xlab = "psrf (Gamma)")
hist(psrf.omega, xlab = "psrf (Omega)")

# Others hist
par(mfrow=c(1,4))

# Beta parameters (species-environment) 
hist(ess.beta, main="ess(beta)")
hist(psrf.beta, main="psrf(beta)")
vioplot(psrf.beta,main="psrf(beta)")
vioplot(ess.beta,main="ess(beta)")

# Gamma parameters (trait-environment)
hist(ess.gamma, main="ess(gamma)")
hist(psrf.gamma, main="psrf(gamma)")
vioplot(psrf.gamma,main="psrf(gamma)")
vioplot(ess.gamma,main="ess(gamma)")

# Omega parameters
hist(ess.omega, main="ess(omega)")
hist(psrf.omega, main="psrf(omega)")
vioplot(psrf.omega,main="psrf(omega)")
vioplot(ess.omega,main="ess(omega)")

## Plot psrf-covariates ----
# Load required packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(stringr)

# Rename columns for easier access
psrf_df <- as.data.frame(psrf.beta)
psrf_df$parameter <- rownames(psrf_df)
rownames(psrf_df) <- NULL
colnames(psrf_df)

colnames(psrf_df) <- c("psrf", "upper_ci", "parameter")

# Extract species names every 10 lines and propagate
species_vector <- psrf_df$parameter[seq(1, nrow(psrf_df), by = 10)] %>%
  str_extract("[^,]+(?=\\]$)")  # everything before the final closing bracket

psrf_data <- psrf_df %>%
  mutate(
    species = rep(
      psrf_df$parameter[seq(1, nrow(psrf_df), by = 10)] %>%
        str_extract("[^,]+(?=\\]$)"),
      each = 10
    ),
    
    # Cleanly extract the variable name
    variable = case_when(
      str_detect(parameter, "poly\\(") ~ paste0(
        str_extract(parameter, "(?<=poly\\().+?(?=,)") %>%
          str_replace_all("_", " "),
        str_extract(parameter, "\\)\\d+") %>% str_remove("\\)")
      ),
      TRUE ~ str_extract(parameter, "(?<=B\\[).+?(?=,)") %>%
        str_replace_all("_", " ") %>%
        str_replace_all("season", "Season ") %>%
        str_trim()
    )
  )

# Flag PSRF > 1.1
psrf_data <- psrf_data %>%
  mutate(flag = psrf > 1.1)

# Create the plot
ggplot(psrf_data, aes(x = variable, y = psrf)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  geom_text_repel(aes(label = ifelse(flag, species, "")), size = 3) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
  labs(
    x = NULL,
    y = "PSRF",
    title = "Convergence of PSRFs by covariate and species"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
