################################################################################
##### Script 1 : Descriptive figures -------------------------------------------
################################################################################

# remove(list=ls())
# see your current working directory
getwd()

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

# We read in the data, set the factors to factors, and have a look at it

fish_survey <- readRDS("./data/fish_survey_Baltic_sea_v2.rds")
fish_survey$Year <- as.factor(fish_survey$Year)
fish_survey$season <- as.factor(fish_survey$season)
data <- fish_survey
head(data)

# Load libraries and functions

library(ggplot2)
library(dplyr)
library(viridis)
library(gridExtra)
library(broom)
library(rlang)
library(mgcv)
library(tidyr)

### Sampling map ----

my_studyArea <- data.frame(lon = data$lon, lat = data$lat)

library("rnaturalearth")
world <- ne_countries(scale = "medium", returnclass = "sf")

q1 <- ggplot(data = world) +
  geom_sf(fill = "grey80", color = "grey10") +
  coord_sf(xlim = c(-12,30), ylim = c(35,71), expand = FALSE)+
  geom_rect(aes(xmin = min(my_studyArea$lon) - 1, xmax = max(my_studyArea$lon) + 1,
                ymin = min(my_studyArea$lat) - 1, ymax = max(my_studyArea$lat) + 1), color = "#D82230", fill = "transparent", linewidth = 1) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white"), # "aliceblue"
        plot.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
q1

q2 <- ggplot(data = world) +
  geom_point(data = my_studyArea, aes(x = lon, y = lat, fill = ""), size = 0.5, color = "#3B5C86") + 
  geom_sf(fill = "grey80", color = "grey10") +
  coord_sf(xlim = range(my_studyArea$lon) + c(-1,1.5), ylim = range(my_studyArea$lat) + c(-1, 1.5), expand = FALSE) +
  guides(fill = guide_legend(title = "Haul", title.position = "right", override.aes = list(size = 4))) +
  labs(x = "", y = "") +
  guides(shape = guide_legend(override.aes = list(size = 0.2)))+
  annotate("text", x = 17.5, y = 56.8, label = "BP", hjust = 0, vjust = 1, size = 8, color = "red") +
  annotate("text", x = 9.8, y = 58.5, label = "S", hjust = 0, vjust = 1, size = 8, color = "red") +
  annotate("text", x = 11, y = 57.1, label = "K", hjust = 0, vjust = 1, size = 8, color = "red") +
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "white"), # "aliceblue"
        plot.margin = unit(c(11, 11, 5.5, 5.5), "pt"),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 8), legend.position = "right")
q2  

library(cowplot)
gg_inset_map2 <-  ggdraw() +
  draw_plot(q2) +
  draw_plot(q1, x = 0.55, y = 0.6, width = 0.5, height = 0.4)

gg_inset_map2

################################################################################
##### Distribution map of a species --------------------------------------------
################################################################################

species_name <- "Clupea harengus"
data_filtered <- data %>% filter(get(species_name) != 0) %>% select(lon, lat)

ggplot(data_filtered, aes(x = lon,y = lat)) + 
  geom_point(size = 0.5, color = "red") +
  borders(fill="grey55", colour = "grey5") +
  theme(
    panel.background = element_rect(fill='lightsteelblue1', colour = 'black'),
    panel.grid.major = element_blank(),) +
  coord_quickmap(xlim=c(min(data$lon), max(data$lon)),ylim= c(min(data$lat),max(data$lat)))+
  labs(title = paste("Distribution map of", species_name), x = "Longitude", y = "Latitude") +
  theme(legend.position = "none")+
  theme(legend.title = element_text(size = 11),legend.text = element_text(size = 10))+
  guides(shape = guide_legend(override.aes = list(size = 1)))

################################################################################
##### Env gradient map ---------------------------------------------------------
################################################################################

### Depth map ----

# Load data
data<- readRDS("./data/fish_survey_Baltic_sea_v2.rds")

# Colour palette
# colpal <- c("midnightblue", "ghostwhite", "coral4")

colpal = viridis(100, option = "plasma")

# Create the graph
ggplot(data, aes(x = lon, y = lat)) + 
  geom_point(aes(color = depth), size = 0.5) +
  scale_color_gradientn(colors = colpal) + 
  borders(fill = "grey55", colour = "grey5") +
  theme(panel.background = element_rect(fill = 'lightsteelblue1', colour = 'black')) +
  coord_quickmap(xlim = c(min(data$lon), max(data$lon)), ylim = c(min(data$lat), max(data$lat))) +
  labs(title = "Depth map", x = "Longitude", y = "Latitude", color = "Depth (m)") +
  theme(legend.position = c(0.4, 0.8),
        legend.background = element_rect(fill = alpha('white', 1), color = 'black')) +
  theme(legend.title = element_text(size = 11), legend.text = element_text(size = 10)) +
  guides(shape = guide_legend(override.aes = list(size = 1)))

ggsave("./plots/Depth map.png", width = 8, height = 6)

### Temperature floor (or other parameter) map ----

data <- fish_survey

# We select to use data from the year 2021 only even if the original data would include repeated surveys to the same locations.
data = droplevels(subset(data,Year==2021))

# We select to use data from the 3rd month only even if the original data would include repeated surveys to the same locations.
data = droplevels(subset(data,Month==3))

col = viridis(100, option = "plasma")

ggplot(data, aes(x = lon, y = lat)) + 
  geom_point(aes(color = T_floor), size = 2) + 
  scale_color_gradientn(colors = col) + 
  borders(fill = "grey55", colour = "grey5") +
  theme(panel.background = element_rect(fill = 'lightsteelblue1', colour = 'black')) +
  coord_quickmap(xlim = c(min(data$lon), max(data$lon)), ylim = c(min(data$lat), max(data$lat))) +
  labs(title = "Bottom temperature map",x = "Longitude", y = "Latitude",  color = "Bottom temperature (°C)") +
  theme(legend.position = c(0.4, 0.75),
        legend.background = element_rect(fill = alpha('white', 1), color = 'black')) +
  theme(legend.title = element_text(size = 11), legend.text = element_text(size = 10)) +
  guides(shape = guide_legend(override.aes = list(size = 1))) 

### Copernicus data depth map ----

depth <- readRDS("./data/Copernicus_BS_extracted/extracted_copernicus_BS_deptho.rds")

ggplot(depth, aes(x = lon, y = lat, fill = depth)) +
  geom_raster() +
  scale_fill_viridis_c() +  
  coord_equal() +
  labs(title = "Depth map", x = "Longitude", y = "Latitude", fill = "Depth (m)") +
  theme_minimal()

### Copernicus data map ----

env <- readRDS("./data/Copernicus_BS_extracted/extracted_copernicus_BS_bottomT.rds")

ggplot(env, aes(x = lon, y = lat, fill = X2021.03.01)) +
  geom_raster() +
  scale_fill_viridis_c() +  
  coord_equal() +
  labs(title = "Bottom temperature map in 03/2021", x = "Longitude", y = "Latitude", fill = "Bottom temperature (°C)") +
  theme_minimal()

### Copernicus map clean----

library(ggplot2)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

env <- readRDS("./data/grid_Copernicus_Baltic_1993-2021/grid_Copernicus_Baltic_2021.rds")
world <- ne_countries(scale = "medium", returnclass = "sf")

lon_limits <- range(env$lon, na.rm = TRUE)
lat_limits <- range(env$lat, na.rm = TRUE)

depth <- ggplot() +
  geom_tile(data = env, aes(x = lon, y = lat, fill = T_floor)) + 
  scale_fill_viridis_c(option = "D") + 
  #scale_fill_gradientn(colours = rev(brewer.pal(11, "RdYlBu"))) +
  geom_sf(data = world, fill = "grey", color = "black", size = 0.2) +
  coord_sf(
    xlim = lon_limits, 
    ylim = lat_limits
  ) +
  labs(
    title = "Bottom temperature map in 2021",
    x = "Longitude",
    y = "Latitude",
    fill = "Temperature (°C)"
  ) +
  annotate("text", x = min(lon_limits), y = max(lat_limits), label = "I", size = 8, fontface = "bold", color = "black") + 
  theme_minimal() +
  theme(
    legend.position = c(0.7, 0.8), 
    legend.background = element_rect(fill = "white", color = "black"), 
    #legend.title = element_text(size = 10, face = "bold"), 
    #legend.key.size = unit(0.8, "cm") 
  )

################################################################################
##### Correlation matrix -------------------------------------------------------
################################################################################

### Correlation matrix for the environmental parameters ----

data <- readRDS("./data/fish_survey_Baltic_sea_v2.rds")
data$season <- as.numeric(as.factor(data$season))
cor.all <- as.data.frame(cor(data[,c( 
  #"Sal_surf", "T_surf", "O2_surf", 
  "Sal_floor", "T_floor", "O2_floor", 
  "depth", "chl"
  #"ph", "season"
)]))

Over.7 <- round(cor.all,2)
Over.7[abs(Over.7) < 0.7] <- NA
Over.7

test <- Over.7

toPlot <- round(cor.all,2)
#toPlot[abs(toPlot) < 0.7] <- 0
#plotOrder <-  corrMatOrder(toPlot, order="AOE") # Euclidian distances
#corrplot(as.matrix(toPlot[plotOrder,plotOrder]), method = "color", col = colorRampPalette(c("deepskyblue2",rep("white", 7),"firebrick1"))(20), mar = c(0,0,0,0), type = "lower", tl.col="black", tl.cex = 1)

all(names(toPlot) == rownames(toPlot))
names(toPlot)

library(ggcorrplot)
ggcorrplot(toPlot, show.diag = T, hc.order = F, type = "lower", legend.title = "Correlation",
           outline.col = "lightgrey", colors = c("deepskyblue2", "white","firebrick1"), lab = TRUE) +
  scale_fill_stepsn(name = "Correlation", breaks = c(-1,-0.7,0.7,1), limits = c(-1,1), show.limits = T,
                    colours = c("deepskyblue2","white","firebrick1"))

ggcorrplot(toPlot, show.diag = TRUE, hc.order = FALSE, type = "lower", 
           outline.col = "lightgrey", colors = c("deepskyblue2", "white","firebrick1"), lab = TRUE) +
  scale_fill_gradient2(name = "Correlation", low = "deepskyblue2", mid = "white", high = "firebrick1", 
                       midpoint = 0, limits = c(-1, 1), space = "Lab")

ggcorrplot(toPlot, show.diag = T, hc.order = F, type = "lower", legend.title = "Correlation",
           outline.col = "lightgrey", colors = c("deepskyblue2", "white","firebrick1"), lab = TRUE, digits = 2) +
  scale_fill_stepsn(name = "Correlation", breaks = c(-1,-0.7,0.7,1), limits = c(-1,1), show.limits = T,
                    colours = c("deepskyblue2","white","firebrick1")) +
  theme(axis.text.x = element_text(colour = c("#009E73", "#009E73", "#009E73",
                                              "#E69F00", "#E69F00", "#E69F00",
                                              "#000000", "#000000", "#000000", "#000000")),
        axis.text.y = element_text(colour = c("#009E73", "#009E73", "#009E73",
                                              "#E69F00", "#E69F00", "#E69F00",
                                              "#000000", "#000000", "#000000", "#000000")),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))

# Limit set at 0.7:

cor.check <- as.data.frame(cor(data[,c(
  "T_surf", "T_floor",
  "Sal_surf", "O2_surf", 
  "Sal_floor", "season"
)]))
Over72 <- round(cor.check,2)
Over72[abs(Over72) < 0.7] <- 0
Over72

ggcorrplot(Over72, show.diag = T, hc.order = TRUE, type = "lower", legend.title = "Correlation",
           outline.col = "lightgrey", colors = c("deepskyblue2", "white","firebrick1"), lab = TRUE) +
  scale_fill_stepsn(name = "Correlation", breaks = c(-1,-0.7,0.7,1), limits = c(-1,1), show.limits = T,
                    colours = c("deepskyblue2","white","firebrick1"))

### Correlation matrix for the fish traits ----

fish_traits <- readRDS("./data/fish_traits.rds")

# Transformer les traits factoriels en variables numériques
alltraits_numeric <- fish_traits
alltraits_numeric$habitat <- as.numeric(as.factor(alltraits_numeric$habitat))
alltraits_numeric$body.shape <- as.numeric(as.factor(alltraits_numeric$body.shape))
alltraits_numeric$feeding.mode <- as.numeric(as.factor(alltraits_numeric$feeding.mode))
alltraits_numeric$fin.shape <- as.numeric(as.factor(alltraits_numeric$fin.shape))
alltraits_numeric$spawning.type <- as.numeric(as.factor(alltraits_numeric$spawning.type))

cor.tr <- as.data.frame(cor(alltraits_numeric[,c( 
  "habitat", "body.shape", "feeding.mode", "fin.shape", "spawning.type",
  "log_age.maturity", "log_length.max", "log_offspring.size",
  "log_fecundity", "log_growth"
)]))

tr.7 <- round(cor.tr,2)
tr.7[abs(tr.7) < 0.7] <- 0
tr.7

library(ggcorrplot)
ggcorrplot(tr.7, show.diag = T, hc.order = TRUE, type = "lower", legend.title = "Correlation",
           outline.col = "lightgrey", colors = c("deepskyblue2", "white","firebrick1"), lab = TRUE) +
  scale_fill_stepsn(name = "Correlation", breaks = c(-1,-0.7,0.7,1), limits = c(-1,1), show.limits = T,
                    colours = c("deepskyblue2","white","firebrick1"))

test <- tr.7
toPlot <- round(cor.tr,2)
all(names(toPlot) == rownames(toPlot))
names(toPlot)

ggcorrplot(toPlot, show.diag = T, hc.order = F, type = "lower", legend.title = "Correlation",
           outline.col = "lightgrey", colors = c("deepskyblue2", "white","firebrick1"), lab = TRUE) +
  scale_fill_stepsn(name = "Correlation", breaks = c(-1,-0.7,0.7,1), limits = c(-1,1), show.limits = T,
                    colours = c("deepskyblue2","white","firebrick1"))

ggcorrplot(toPlot, show.diag = TRUE, hc.order = FALSE, type = "lower", 
           outline.col = "lightgrey", colors = c("deepskyblue2", "white","firebrick1"), lab = TRUE) +
  scale_fill_gradient2(name = "Correlation", low = "deepskyblue2", mid = "white", high = "firebrick1", 
                       midpoint = 0, limits = c(-1, 1), space = "Lab")

