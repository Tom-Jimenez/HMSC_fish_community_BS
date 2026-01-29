################################################################################
##### Script 5 : Explore parameter estimates -----------------------------------
################################################################################

# remove(list=ls())

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

################################################################################
##### Variance partitioning ----------------------------------------------------
################################################################################

head(m$X)

library(viridis)
par(mfrow=c(1,1))
groupnames = c("T_floor", "Sal_floor", "depth", "O2_floor", "chl", "season")
               
group = c(1,1,2,2,3,4,5,6,6)
VP = computeVariancePartitioning(m, group = group, groupnames = groupnames)

library(RColorBrewer)  
cols <- brewer.pal(8, "Paired")

#windows(3000, 2400, pointsize = 12)
plotVariancePartitioning(m,VP, cols=cols)

par(mar=c(9,4,1,1)+.1)
plotVariancePartitioning(m,VP, cex.names = 0.7, las=2, mar = c(5, 8, 4, 2), viridis(8))

plotVariancePartitioning(m,VP, cex.names = 0.7, las=2, mar = c(5, 8, 4, 2), viridis(8),legend.text = FALSE)

dev.off()

plotVariancePartitioning <- function (hM, VP, cols = NULL, main = "", 
          ...) 
{
  ng = dim(VP$vals)[1]
  if (is.null(cols)) {
    cols = heat.colors(ng, alpha = 1)
  }
  leg = VP$groupnames
  for (r in 1:hM$nr) {
    leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
  }
  means = round(100 * rowMeans(VP$vals), 1)
  for (i in 1:ng) {
    leg[i] = paste(leg[i], " (mean = ", toString(means[i]), 
                   ")", sep = "")
  }
  barplot(VP$vals, main = main, xlab = "", ylab = "Variance proportion", 
          las = 1, legend = leg, col = cols, ...)
}

par(mar = c(14,4,1,1) + .1)
plotVariancePartitioning(
  m, VP,
  cols = cols,
  cex.names = 1.2,
  las = 2,
  mar = c(5, 8, 4, 2),
  legend.text = F
)

dev.off()

################################################################################
##### Beta plot ----------------------------------------------------------------
################################################################################

# We next construct a beta-plot showing the estimates of species niche parameters

postBeta = getPostEstimate(m, parName="Beta")
plotBeta(m, post=postBeta, plotTree = TRUE, spNamesNumbers=c(FALSE, FALSE))

##### Beta plot v2 -------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggtree)
library(ggimage)
library(patchwork)

postBeta = getPostEstimate(m, parName="Beta")
supportLevel <- 0.95
plotBeta(m, postBeta, main = paste0("Beta\n", "support level ", supportLevel), 
         supportLevel = supportLevel, param = "Sign", 
         plotTree = FALSE, spNamesNumbers = c(T,F), covNamesNumbers = c(T, FALSE), 
         colors = colorRampPalette(c("#003366","white","#E31B23")))

load("./data/tree_plot.RData")
tree_plot$tip.label <- gsub("_", " ", tree_plot$tip.label)
setdiff(colnames(m$Y), tree_plot$tip.label)

tree_plot <- ggtree(tree_plot, branch.length="none") + geom_tiplab(fontface = "italic") + xlim(0, 30)

class(tree_plot)
length(tree_plot$edge.length)
nrow(tree_plot$edge)
# un objet phylo doit avoir autant de longueurs que d’arêtes mais ce n'est pas le cas ici

library(ape)
tree_plot <- reorder.phylo(tree_plot, "postorder")
ggtree(tree_plot, branch.length = "none") +
  geom_tiplab(fontface = "italic") +
  xlim(0, 30)
tree_plot

postBeta <-  getPostEstimate(m, parName = "Beta")

temp.beta <- as.data.frame(postBeta$mean) # 'temp' stands for temporary
temp.beta$param <- m$covNames
temp.beta 
temp.beta <- temp.beta %>% tidyr::pivot_longer(.,-param, names_to = 'species', values_to ='value')

sup.beta <- as.data.frame(postBeta$support) # n.beta because I'm extracting the negative values
sup.beta$param <- m$covNames
sup.beta 
sup.beta <- sup.beta %>% tidyr::pivot_longer(.,-param, names_to = 'species', values_to ='value')

temp.beta$support <- sup.beta$value

supportLevel <- 0.95
temp.beta$pos.neg <- 0
temp.beta$pos.neg[temp.beta$support > supportLevel] <- 1
temp.beta$pos.neg[temp.beta$support < (1-supportLevel)] <- -1

temp.beta$pos.neg <- factor(temp.beta$pos.neg, levels = c(1,0,-1)) # this is needed to always plot the '+' and '-' in the legend

temp.beta$Environment <- factor(temp.beta$param)
levels(temp.beta$Environment)
temp.beta$Environment <- dplyr::recode(temp.beta$Environment, 
                                       "chl" = "chl",
                                       "depth" = "depth",
                                       "poly(T_floor, degree = 2, raw = TRUE)1" = "T_floor1",
                                       "poly(T_floor, degree = 2, raw = TRUE)2" = "T_floor2",
                                       "poly(Sal_floor, degree = 2, raw = TRUE)1" = "Sal_floor1",
                                       "poly(Sal_floor, degree = 2, raw = TRUE)2" = "Sal_floor2",
                                       "O2_floor" = "O2_floor",
                                       "seasonSummer" = "season Summer",
                                       "seasonWinter" = "season Winter")
temp.beta$Environment <-  factor(temp.beta$Environment,
                                 levels = c("(Intercept)", "depth", "Sal_floor1", "Sal_floor2",
                                            "T_floor1", "T_floor2", "O2_floor", "chl",
                                            "season Summer", "season Winter"))
levels(temp.beta$Environment)

# Set species order same as phylogenetic tree

sp_order <- c(
  "Glyptocephalus cynoglossus",
  "Microstomus kitt",
  "Platichthys flesus",
  "Limanda limanda",
  "Pleuronectes platessa",
  "Hippoglossoides platessoides",
  "Solea solea",
  "Buglossidium luteum",
  "Scophthalmus rhombus",
  "Arnoglossus laterna",
  "Hyperoplus immaculatus",
  "Hyperoplus lanceolatus",
  "Ammodytes marinus",
  "Chelidonichthys lucerna",
  "Eutrigla gurnardus",
  "Myoxocephalus scorpius",
  "Trachinus draco",
  "Agonus cataphractus",
  "Lumpenus lampretaeformis",
  "Pollachius pollachius",
  "Pollachius virens",
  "Trisopterus minutus",
  "Trisopterus esmarkii",
  "Melanogrammus aeglefinus",
  "Gadus morhua",
  "Merlangius merlangus",
  "Merluccius merluccius",
  "Enchelyopus cimbrius",
  "Sprattus sprattus",
  "Clupea harengus",
  "Engraulis encrasicolus",
  "Myoxocephalus quadricornis",
  "Gasterosteus aculeatus",
  "Aphia minuta",
  "Mullus surmuletus",
  "Osmerus eperlanus",
  "Trachurus trachurus",
  "Entelurus aequoreus",
  "Scomber scombrus",
  "Callionymus lyra"
)

temp.beta$Species_sort <- as.factor(temp.beta$species)
temp.beta$Species_sort <- factor(temp.beta$species, levels = rev(sp_order))
temp.beta <- temp.beta[order(temp.beta$Species_sort), ]
levels(temp.beta$Species_sort)

beta.heatmap <- ggplot(temp.beta, aes(x = Environment, y = Species_sort)) +
  geom_tile(aes(fill=pos.neg),color='black') +
  scale_fill_manual(values=c('#E31B23','white','#003366'),labels=c('+','','-'), name='',guide=guide_legend(keyheight =4, keywidth = 1)) + 
  xlab('') + ylab('')+
  theme_bw()+
  theme(text = element_text(family = "sans",size = 12),
        strip.text = element_text(family = "sans",face = 'bold'),
        axis.text = element_text(family = "sans"),
        legend.title = element_text(family = "sans"),
        legend.text = element_text(family = "sans"),axis.text.x  = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.grid = element_blank(),
        panel.background = element_blank(), strip.background = element_blank())

beta.heatmap_nonames <- ggplot(temp.beta, aes(x = Environment, y = Species_sort)) +
  geom_tile(aes(fill=pos.neg),color='black')+
  scale_fill_manual(values=c('#E31B23','white','#003366'),labels = c('+','','-'), name = '', guide = guide_legend(keyheight = 4, keywidth = 1)) +
  theme_bw()+
  labs(y = "", x = "") +
  theme(text = element_text(family = "sans", size = 10),
        strip.text = element_text(family = "sans",face = 'bold'),
        axis.text = element_text(family = "sans"),
        legend.title = element_text(family = "sans"),
        legend.text = element_text(family = "sans"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.grid = element_blank(),
        panel.background = element_blank(), strip.background = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_blank())

# Check if species names match
tree_plot + beta.heatmap + plot_layout(widths = c(2, 3))
tree_plot$data$label <- paste0(" ", tree_plot$data$label) # add space between tree and sp names
tree_plot + beta.heatmap + plot_layout(widths = c(2, 3))

tree_plot + beta.heatmap_nonames + plot_layout(widths = c(3, 4))

tree_plot <- ggtree(tree_plot, branch.length = "none") + geom_tiplab(fontface = "italic", size = 3) + xlim(0, 30)

# If they match, continue
#windows(3000, 2400, pointsize = 12)
tree_plot + beta.heatmap + plot_layout(widths = c(4.5, 5))
tree_plot + beta.heatmap_nonames + plot_layout(widths = c(4.5, 5))

# Add proportion of responses at the bottom
betaSummary <- ggplot(temp.beta, aes(x = Environment, fill = pos.neg)) + 
  geom_bar(colour = "black") +
  scale_fill_manual(values = c('#E31B23','white','#003366'), labels = c('+','','-'),
                    name = '',guide = guide_legend(keyheight = 4, keywidth = 1)) +
  labs(y = "Species", x = "") +
  theme_bw() +
  theme(text = element_text(family = "sans",size = 12),
        strip.text = element_text(family = "sans",face = 'bold'),
        axis.text = element_text(family = "sans"),
        legend.position = "none",
        legend.title = element_text(family = "sans"),
        legend.text = element_text(family = "sans"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.grid = element_blank(),
        panel.background = element_blank(), strip.background = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(family = "sans"),
        plot.margin = unit(c(0,1,1,1), "cm"))
betaSummary

dev.off()

################################################################################
##### Gamma plot ---------------------------------------------------------------
################################################################################

# We examine next if the species niches are linked to their traits with a Gamma-plot

postGamma = getPostEstimate(m, parName="Gamma")
supportLevel <- 0.95
plotGamma(m, postGamma, main = paste0("Gamma\n", "support level ", supportLevel), 
          supportLevel = supportLevel, param = "Sign", covNamesNumbers = c(T, FALSE), 
          colors = colorRampPalette(c("#003366","white","#E31B23")))

# Let us next relax the level of required statistical support from 0.95 to 0.85.

supportLevel <- 0.85
plotGamma(m, postGamma, main = paste0("Gamma\n", "support level ", supportLevel), 
          supportLevel = supportLevel, param = "Sign", covNamesNumbers = c(T, FALSE), 
          colors = colorRampPalette(c("#003366","white","#E31B23")))

##### Gamma plot v2 ------------------------------------------------------------

library(ggimage)
library(ggtree)
library(dplyr)

supportLevel <- 0.85

# gamma heatmap
postGamma <-  getPostEstimate(m, parName = "Gamma")
temp.gamma <- as.data.frame(postGamma$mean) # 'temp' stands for temporary
names(temp.gamma) <- m$trNames
temp.gamma$param <- m$covNames
temp.gamma 
temp.gamma <- temp.gamma %>% tidyr::pivot_longer(.,-param, names_to = 'traits', values_to ='value')

sup.gamma <- as.data.frame(postGamma$support) # n.beta because I'm extracting the negative values
names(sup.gamma) <- m$trNames
sup.gamma$param <- m$covNames
sup.gamma 
sup.gamma <- sup.gamma %>% tidyr::pivot_longer(.,-param, names_to = 'traits', values_to ='value')

temp.gamma$support <- sup.gamma$value

temp.gamma$pos.neg <- 0
temp.gamma$pos.neg[temp.gamma$support > supportLevel] <- 1
temp.gamma$pos.neg[temp.gamma$support < (1-supportLevel)] <- -1

temp.gamma$pos.neg <- factor(temp.gamma$pos.neg, levels = c(1,0,-1)) # this is needed to always plot the '+' and '-' in the legend

# Rename variables
temp.gamma$Environment <- factor(temp.gamma$param)
levels(temp.gamma$Environment)
temp.gamma$Environment <- dplyr::recode(temp.gamma$Environment, 
                                        "chl" = "chl",
                                        "depth" = "depth",
                                        "poly(T_floor, degree = 2, raw = TRUE)1" = "T_floor1",
                                        "poly(T_floor, degree = 2, raw = TRUE)2" = "T_floor2",
                                        "poly(Sal_floor, degree = 2, raw = TRUE)1" = "Sal_floor1",
                                        "poly(Sal_floor, degree = 2, raw = TRUE)2" = "Sal_floor2",
                                        "O2_floor" = "O2_floor",
                                        "seasonSummer" = "season Summer",
                                        "seasonWinter" = "season Winter")
temp.gamma$Environment <-  factor(temp.gamma$Environment,
                                  levels = c("(Intercept)", "depth", "Sal_floor1", "Sal_floor2",
                                             "T_floor1", "T_floor2", "O2_floor", "chl",
                                             "season Summer", "season Winter"))
levels(temp.gamma$Environment)

#windows(800, 600, pointsize = 12) #opens a separate window with the size you want
par(mar = c(4, 12, 0, 0))
plotGamma(m, postGamma, main = paste0("Gamma\n", "support level ", supportLevel), 
          supportLevel = supportLevel, param = "Sign", covNamesNumbers = c(T, FALSE), 
          colors = colorRampPalette(c("#003366","white","#E31B23")))

# Rename traits:
temp.gamma$Traits <- factor(temp.gamma$traits)
levels(temp.gamma$Traits)
temp.gamma$Traits <- dplyr::recode(temp.gamma$Traits, 
                                   "body.shapeelongated" = "Elongated (body shape)",
                                   "body.shapeflat" = "Flat  (body shape)",
                                   "body.shapefusiform" = "Fusiform  (body shape)",
                                   "feeding.modegeneralist/piscivorous" =  "Generalist / Piscivorous (feeding mode)",
                                   "feeding.modeplanktivorous" = "Planktivorous (feeding mode)",
                                   "habitatdemersal" = "Demersal (habitat)",
                                   "habitatpelagic" = "Pelagic (habitat)",
                                   "log_age.maturity" = "Age at maturity",
                                   "log_fecundity" = "Fecundity",
                                   "log_growth" = "Growth (K)",
                                   "log_length.max" = "Maximum length",
                                   "log_offspring.size" = "Offspring size"
)
sort(levels(temp.gamma$Traits))
temp.gamma$Traits <-  factor(temp.gamma$Traits,
                             levels = rev(c("Age at maturity", "Fecundity", "Growth (K)", "Maximum length", "Offspring size",
                                            "Elongated (body shape)", "Flat  (body shape)", "Fusiform  (body shape)",
                                            "Generalist / Piscivorous (feeding mode)", "Planktivorous (feeding mode)",
                                            "Demersal (habitat)", "Pelagic (habitat)",
                                            "(Intercept)"
                             )))
levels(temp.gamma$Traits)

#windows(800, 600, pointsize = 12) #opens a separate window with the size you want
ggplot(temp.gamma, aes(x = Environment, y = Traits)) +
  geom_tile(aes(fill=pos.neg),color='black')+
  scale_fill_manual(values=c('#E31B23','white','#003366'),labels=c('+','','-'), 
                    name='',guide=guide_legend(keyheight =4, keywidth = 1)) + 
  xlab('') + ylab('')+
  theme_bw()+
  theme(text = element_text(family = "sans",size = 14),
        strip.text = element_text(family = "sans",face = 'bold'),
        axis.text = element_text(family = "sans"),
        legend.title = element_text(family = "sans"),
        legend.text = element_text(family = "sans"),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.grid = element_blank(),
        panel.background = element_blank(), strip.background = element_blank())

# Add grouping 

gamma.heatmap <- ggplot(temp.gamma, aes(x = Environment, y = Traits)) +
  geom_tile(aes(fill=pos.neg),color='black')+
  scale_fill_manual(values=c('#E31B23','white','#003366'),labels=c('+','','-'), name='',guide=guide_legend(keyheight = 4, keywidth = 1)) + xlab('') + ylab('')+
  theme_bw()+
  theme(text = element_text(family = "sans",size = 10),
        strip.text = element_text(family = "sans",face = 'bold'),
        axis.text = element_text(family = "sans"),
        legend.title = element_text(family = "sans"),
        legend.text = element_text(family = "sans"),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y  = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.grid = element_blank(),
        panel.background = element_blank(), strip.background = element_blank())
gamma.heatmap
dev.off()

myGamma <- temp.gamma
myGamma$Group <- stringr::str_extract(string = myGamma$Traits, pattern = "(?<=\\().*(?=\\))")
myGamma$Group[myGamma$Group == "K"] <- NA

myGamma$Only_traits <- gsub("\\s*\\([^\\)]+\\)", "", myGamma$Traits)
myGamma$Only_traits[myGamma$Only_traits == "Growth"] <- "Growth (K)"
myGamma$Only_traits[myGamma$Group == "Intercept"] <- "Intercept"
#myGamma$Group[myGamma$Group == "Intercept"] <- ""
myGamma$Group[is.na(myGamma$Group)] <- ""
myGamma$Group <- gsub(" ", "\n", myGamma$Group)
myGamma$Group <- stringr::str_to_title(myGamma$Group) 
myGamma$Group[myGamma$Group == "Caudal\nFin\nShape"] <- "Caudal Fin\nShape"

# Add grouping to continuous variables:
unique(myGamma$Group)
myGamma$Group[myGamma$Group == ""] <- "Continuous\nTraits"

# Set level order
myGamma$Group <- factor(myGamma$Group, levels = c("Intercept", "Continuous\nTraits", "Body\nShape", "Caudal Fin\nShape", "Feeding\nMode", "Habitat", "Spawning\nType"))
myGamma$Only_traits <- factor(myGamma$Only_traits, levels = rev(sort(unique(myGamma$Only_traits))))

gamma.plot.opt2 <- ggplot(myGamma, aes(x = Environment, y = Only_traits, fill = pos.neg)) + 
  geom_tile(color = 'black')+
  scale_fill_manual(values = c('#E31B23','white','#003366'), labels = c('+','','-'), 
                    name='',guide = guide_legend(keyheight = 4, keywidth = 1)) +
  facet_grid(Group ~., scales = "free_y", space = "free_y") +
  xlab('') + ylab('') +
  theme_bw() +
  theme(text = element_text(family = "sans", size = 14),
        axis.text = element_text(family = "sans"),
        legend.title = element_text(family = "sans"),
        legend.text = element_text(family = "sans", size = 16),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y  = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(family = "sans",face = 'bold'),
        strip.text.y = element_text(size = 12, face = "bold", hjust = 0.5, angle = 0),
        strip.placement = "outside", strip.background = element_blank())
gamma.plot.opt2

dev.off()

# Another way of examing the influence of traits is to see how much of the variation they
# explain among the responses of the species to their covariates.
VP$R2T$Beta
VP$R2T$Y

################################################################################
##### Species associations -----------------------------------------------------
################################################################################

# We next evaluate the posterior distribution of the phylogenetic signal in species niches

mpost = convertToCodaObject(m)
round(summary(mpost$Rho, quantiles = c(0.025, 0.5, 0.975))[[2]],2)

# We next illustrate the species associations revealed by the random effects with the corrplot function.

library(ggcorrplot)
library(corrplot)
OmegaCor = computeAssociations(m)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
toPlot = sign(toPlot)
plotOrder = corrMatOrder(OmegaCor[[1]]$mean,order="AOE")
corrplot(toPlot[plotOrder,plotOrder], method = "color", tl.cex=0.5,
         col=colorRampPalette(c("blue", "white", "red"))(255))

# Species associations

supportLevel <- 0.95
OmegaCor <- computeAssociations(m)
names(OmegaCor) <- names(m$ranLevels)
toPlot <- ((OmegaCor$cell$support > supportLevel) + 
             (OmegaCor$cell$support < (1 - supportLevel)) > 0) * OmegaCor$cell$mean 
plotOrder <- corrMatOrder(OmegaCor$cell$mean, order="AOE")

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

corr <- as.matrix(toPlot[plotOrder,plotOrder])
# corr <- base::round(x = corr, digits = 2)
corr[corr > 0] <- 1
corr[corr < 0] <- -1
upper_tri <- get_upper_tri(corr)

melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)
melted_cormat$value <- factor(melted_cormat$value, levels = c("1", "0", "-1"))

# Heatmap
omega.plot <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) + 
  scale_x_discrete(limits = rev(levels(melted_cormat$Var2)))+
  geom_tile(color = "black") +
  labs(x = "", y = "") +
  scale_fill_manual(values = c('#E31B23','white','#003366'), 
                    labels = c('+','','-'), name ='', 
                    guide = guide_legend(keyheight = 4, keywidth = 1)) +
  theme_bw() +
  theme(text = element_text(family = "sans", size = 10),
        strip.text = element_text(family = "Helvetica", face = 'bold'),
        axis.text = element_text(face = "italic", color = "black", size = 8),
        legend.text = element_text(size = 10),
        axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.25),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.background = element_blank(), strip.background = element_blank(),
        legend.position = c(0.96,0.6),
        plot.margin = unit(c(11, 11, 5.5, 5.5), "pt"))
omega.plot

dev.off()

# We next examine at which spatial scale the variation captured by the random effect occurs.
mpost = convertToCodaObject(m) 
summary(mpost$Alpha[[1]], quantiles = c(0.025, 0.5, 0.975))

################################################################################
##### Proportion of the response to covariates explained by traits -------------
################################################################################

round(head(m$X),2)
colnames(m$X)
groupnames = c("T_floor", "Sal_floor", "depth", "O2_floor", "chl", "season")
group = c(1,1,1,2,2,3,4,5,6,6)

VPfullPa <- computeVariancePartitioning(m, group = group, groupnames = groupnames)
trait.beta <- data.frame(FullPa = c(VPfullPa$R2T$Beta*100))
trait.beta <- round(trait.beta, 2)
trait.beta$Environment <- c("(Intercept)", "depth", "Sal_floor1", "Sal_floor2",
                            "T_floor1", "T_floor2", "O2_floor", "chl", "season Summer", "season Winter")
trait.beta

ggplot(trait.beta, aes(x = reorder(Environment, - FullPa), y = FullPa)) + geom_bar(stat = "identity", position = position_dodge(.9), fill = "skyblue4") +
  labs(subtitle = "Proportion of response to covariates explained by traits", y = "Proportion (%)", x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = c(0.8, 0.6),
        text = element_text(size = 18))

dev.off()

################################################################################
##### Latent factors -----------------------------------------------------------
################################################################################

### Latent variables (cell)

postEta <- getPostEstimate(m, parName = "Eta")
dim(postEta$mean)

postEtamean1 <- postEta$mean[,1] # Latent factor 1
postEtamean2 <- postEta$mean[,2] # Latent factor 2
postEtamean3 <- postEta$mean[,3] # Latent factor 3
postEtamean4 <- postEta$mean[,4] # Latent factor 4
postEtamean5 <- postEta$mean[,5] # Latent factor 5
postEtamean6 <- postEta$mean[,6] # Latent factor 6

xy_cell = data.frame(cell = data$cell, x_coordinate = data$Lon_c, y_coordinate = data$Lat_c)
xy_cell_unique <- xy_cell[!duplicated(xy_cell), ]
dataLf <- cbind(as.data.frame(xy_cell_unique), postEtamean1, postEtamean2, postEtamean3, postEtamean4, postEtamean5, postEtamean6)

range(dataLf$postEtamean1)
range(dataLf$postEtamean2)
range(dataLf$postEtamean3)
range(dataLf$postEtamean4)
range(dataLf$postEtamean5)
range(dataLf$postEtamean6)
lv_range <- range(dataLf[,c("postEtamean1", "postEtamean2", "postEtamean3", "postEtamean4", "postEtamean5", "postEtamean6")])

library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) + geom_sf(fill = "grey80", color = "grey20") +
  coord_sf(xlim = range(dataLf$'x_coordinate') + c(-1,1.5), ylim = range(dataLf$'y_coordinate') + c(-1,1.5), expand = FALSE)

library(paletteer)

latfac1 <- ggplot(data = world) + 
  geom_point(data = dataLf, aes(x = `x_coordinate`, y = `y_coordinate`, color = postEtamean1), size = 3) +
  geom_sf(fill = "grey80", color = "grey20") +
  coord_sf(xlim = range(dataLf$`x_coordinate`) + c(-1, 1), 
           ylim = range(dataLf$`y_coordinate`) + c(-1, 1), 
           expand = FALSE) +
  scale_color_gradient2(low = "#003366", mid = "white", high = "#E31B23", 
                        midpoint = mean(dataLf$postEtamean1), 
                        name = "Latent\nfactor 1") +
  annotate("text", x = min(dataLf$`x_coordinate`), 
           y = max(dataLf$`y_coordinate`) + 0.5, 
           label = "A", size = 6, fontface = "bold", color = "black") + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

latfac2 <- ggplot(data = world) + 
  geom_point(data = dataLf, aes(x = `x_coordinate`, y = `y_coordinate`, color = postEtamean2), size = 3) +
  geom_sf(fill = "grey80", color = "grey20") +
   coord_sf(xlim = range(dataLf$`x_coordinate`) + c(-1, 1), 
           ylim = range(dataLf$`y_coordinate`) + c(-1, 1), 
           expand = FALSE) +
  scale_color_gradient2(low = "#003366", mid = "white", high = "#E31B23", 
                        midpoint = mean(dataLf$postEtamean2), 
                        name = "Latent\nfactor 2") +
  annotate("text", x = min(dataLf$`x_coordinate`), 
           y = max(dataLf$`y_coordinate`) + 0.5, 
           label = "B", size = 6, fontface = "bold", color = "black") + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

latfac3 <- ggplot(data = world) + 
  geom_point(data = dataLf, aes(x = `x_coordinate`, y = `y_coordinate`, color = postEtamean3), size = 3) +
  geom_sf(fill = "grey80", color = "grey20") +
    coord_sf(xlim = range(dataLf$`x_coordinate`) + c(-1, 1), 
           ylim = range(dataLf$`y_coordinate`) + c(-1, 1), 
           expand = FALSE) +
  scale_color_gradient2(low = "#003366", mid = "white", high = "#E31B23", 
                        midpoint = mean(dataLf$postEtamean3), 
                        name = "Latent\nfactor 3") +
  annotate("text", x = min(dataLf$`x_coordinate`), 
           y = max(dataLf$`y_coordinate`) + 0.5, 
           label = "C", size = 6, fontface = "bold", color = "black") + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
