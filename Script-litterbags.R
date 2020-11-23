###################################
# DCE Mangroves
###################################

###################################
# Packages
###################################
library(ade4)
library(factoextra)
library(cowplot)
library(vegan)

###################################
# Data pigments
###################################

percent <- read.table("litterbags-pigments-all.txt",h=T)
PERCENT <- read.table("litterbags-pigments-short.txt",h=T)

###################################
# Aesthetics
###################################

my.palette <- colorRampPalette(c("red3","orange","yellow","green3","royalblue"))

###################################
# ACP pigments
###################################

names(PERCENT)
res.pca <- dudi.pca(PERCENT,scannf=F,nf=3)

res.acp.var <- get_pca_var(res.pca) # obtenir les cos2/corel/contrib des variables
res.acp.ind <- get_pca_ind(res.pca) # obtenir les cos2/contrib des individus

VAR <- fviz_pca_var(res.pca, col.var="cos2", # colorier les variables en fonction de leur contrib
             gradient.cols = my.palette(3), # controler les couleurs
             repel = TRUE # éviter l'overlapping
)

fac <- factor(paste(percent$station))
IND <- fviz_pca_ind(res.pca,
                    label = "none", # cacher les étiquettes
                    habillage = fac, # couleur selon le facteur
                    palette = my.palette(length(levels(fac))), # controler les couleurs
                    addEllipses = TRUE # ajouter des ellipses
)

plot_grid(IND,VAR,labels=c("a","b"),ncol=1)

###################################
# BCA pigments
###################################

res.bca <- bca(res.pca,fac,scannf=F,nf=3)

ggplot(res.bca$ls,
  aes(x=CS1,y=CS2,col=fac))+
  geom_point(size=3)+
  stat_ellipse()+
  theme_bw()

  
  
  
  
  

###################################
# PCoA pigments
###################################

res.pcoa <- dudi.pco(cailliez(vegdist(PERCENT,"bray")),scannf=F,nf=3)



plot_grid(IND,VAR,labels=c("a","b"),ncol=1)

