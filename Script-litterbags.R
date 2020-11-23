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

names(percent)
res.pca <- dudi.pca(PERCENT,scannf=F,nf=3)

res.acp.var <- get_pca_var(res.pca) # obtenir les cos2/corel/contrib des variables
res.acp.ind <- get_pca_ind(res.pca) # obtenir les cos2/contrib des individus

VAR <- fviz_pca_var(res.pca, col.var="cos2", # colorier les variables en fonction de leur contrib
             gradient.cols = my.palette(3), # controler les couleurs
             repel = TRUE # éviter l'overlapping
)

fac <- factor(paste(percent$station))
lab<- factor(paste(percent$time))
IND <- fviz_pca_ind(res.pca,
                    label = "none",
                    habillage = fac, # couleur selon le facteur
                    addEllipses = TRUE # ajouter des ellipses
)

plot_grid(IND,VAR,labels=c("a","b"),ncol=1)

###################################
# BCA pigments
###################################

res.bca <- bca(res.pca,fac,scannf=F,nf=3)
res.bca$ratio
IND.bca <- ggplot(res.bca$ls,
  aes(x=CS1,y=CS2,col=fac))+
  geom_point(size=3)+
  stat_ellipse()+
  theme_linedraw(base_size = 14)

VAR.bca <- ggplot(res.bca$co,aes(x=Comp1,y=Comp2,label=rownames(res.bca$co)))+
  geom_text(col="red")+
  geom_segment(aes(x=0,y=0,xend=Comp1,yend=Comp2),
               arrow = arrow(length = unit(0.5, "cm")),
               col="grey")+
  theme_linedraw(base_size = 14)

plot_grid(IND.bca,VAR.bca,labels=c("a","b"),ncol=1)



###################################
# PCoA pigments
###################################

res.pcoa <- dudi.pco(cailliez(vegdist(PERCENT,"bray")),
                     scannf=F,
                     nf=3)


ggplot(res.pcoa$li,
       aes(x=A1,y=A2,col=fac))+
  geom_point(size=3)+
  stat_ellipse()+
  theme_linedraw(base_size = 14)


###################################
# nMDS pigments
###################################

res.nmds <- metaMDS(PERCENT,"bray",k=2)
summary(res.nmds)

ggplot(data.frame(res.nmds$points),
       aes(x=MDS1,y=MDS2,col=fac))+
  geom_point(size=3)+
  stat_ellipse()+
  theme_linedraw(base_size = 14)





