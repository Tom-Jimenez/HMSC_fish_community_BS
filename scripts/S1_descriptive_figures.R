################################################################################
##### Script 1 : Descriptive figures -------------------------------------------
##### Tom Jimenez, 16-11-2026 --------------------------------------------------
################################################################################

remove(list=ls())
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

# Load libraries

library(ggplot2)
library(dplyr)
library(viridis)
library(gridExtra)
library(broom)
library(rlang)
library(mgcv)
library(tidyr)

################################################################################
##### Sampling map -------------------------------------------------------------
################################################################################

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
#corrplot(as.matrix(toPlot[plotOrder,plotOrder]), method = "color", col = colorRampPalette(c("#2166AC",rep("white", 7),"#E31B23"))(20), mar = c(0,0,0,0), type = "lower", tl.col="black", tl.cex = 1)

all(names(toPlot) == rownames(toPlot))
names(toPlot)

library(ggcorrplot)
ggcorrplot(toPlot, show.diag = T, hc.order = F, type = "lower", legend.title = "Correlation",
           outline.col = "lightgrey", colors = c("#2166AC", "white","#E31B23"), lab = TRUE) +
  scale_fill_stepsn(name = "Correlation", breaks = c(-1,-0.7,0.7,1), limits = c(-1,1), show.limits = T,
                    colours = c("#2166AC","white","#E31B23"))

ggcorrplot(toPlot, show.diag = TRUE, hc.order = FALSE, type = "lower", 
           outline.col = "lightgrey", colors = c("#2166AC", "white","#E31B23"), lab = TRUE) +
  scale_fill_gradient2(name = "Correlation", low = "#2166AC", mid = "white", high = "#E31B23", 
                       midpoint = 0, limits = c(-1, 1), space = "Lab")

dev.off()
