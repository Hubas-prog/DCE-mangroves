###################################
# DCE Mangroves
###################################

###################################
# Packages
###################################
library(ade4)
library(vegan)
library(factoextra)
library(cowplot)
library(viridis)

###################################
# Data pigments + data microbiote
###################################

setwd("~/sync/mangroves/DCE-mangroves/")

pig <- read.table("litterbags-pigments-all.txt",h=T)
pig_short <- read.table("litterbags-pigments-short.txt",h=T)
pig$ID <- paste(substr(pig$ID,3,4), substr(pig$ID,1,2), substr(pig$ID,5,6), sep="")
rownames(pig) <- pig$ID
dim(pig)

mic16S <- read.csv2("table_OTU_table97_guyane_LB_16SV4.csv", header=T,sep="\t")
mic16S <- as.data.frame(t(mic16S))
dim(mic16S)
colnames(mic16S) <- paste0("OTU", 1:ncol(mic16S))
mic16S <- mic16S[-1:-2,]
mic16S[1:10,1:10]

rownames(pig)
rownames(mic16S) -> mic16S$ID

table1 <- merge(pig,mic16S, by="ID")
dim(table1)
table1[1:10,25:35]








###################################
# Aesthetics
###################################

my.palette <- colorRampPalette(c("red3","orange","green3","royalblue"))
fac <- factor(paste(percent$station))
lab<- factor(paste(percent$time))
cross<-factor(paste(fac,lab))

###################################
# plot function
###################################

my_ind_plot<-function(data.li,groups,Xaxis,Yaxis,pal.name){
  plot(data.li[,Yaxis]~data.li[,Xaxis],
       type="n",
       xlab="Axis1",
       ylab="Axis2")
  s.class(data.li[,c(Xaxis,Yaxis)],
          groups,
          add.plot=T,
          axesell=F,
          col=alpha(rep(my.palette(3),each=4),0.7),
          cellipse = 0)
}

###################################
# ACP pigments
###################################

names(percent)
res.pca <- dudi.pca(PERCENT,scannf=F,nf=3)

res.acp.var <- get_pca_var(res.pca) # obtenir les cos2/corel/contrib des variables
res.acp.ind <- get_pca_ind(res.pca) # obtenir les cos2/contrib des individus

fviz_pca_var(res.pca, col.var="cos2", # colorier les variables en fonction de leur contrib
             repel = T# éviter l'overlapping
)

my_ind_plot(res.pca$li,
            groups=cross,
            pal.name="Dark2",
            Xaxis=1,
            Yaxis=2)


###################################
# BCA pigments
###################################

res.bca <- bca(res.pca,fac,scannf=F,nf=3)
res.bca$ratio
IND.bca <- ggplot(res.bca$ls,
  aes(x=CS1,y=CS2,col=fac))+
  geom_point(aes(size=lab))+
  stat_ellipse()+
  theme_linedraw(base_size = 14)

VAR.bca <- ggplot(res.bca$co,aes(x=Comp1,
                                 y=Comp2,
                                 label=rownames(res.bca$co),
                                 ))+
  geom_text(col="red")+
  geom_segment(aes(x=0,y=0,xend=Comp1,yend=Comp2),
               arrow = arrow(length = unit(0.5, "cm")),
               col="grey")+
  theme_linedraw(base_size = 14)

plot_grid(IND.bca,VAR.bca,labels=c("a","b"),ncol=2)



###################################
# PCoA pigments
###################################

res.pcoa <- dudi.pco(cailliez(vegdist(PERCENT,"bray")),
                     scannf=F,
                     nf=3)


ggplot(res.pcoa$li,
       aes(x=A1,y=A2,col=fac))+
  geom_point(aes(size=lab))+
  stat_ellipse()+
  theme_linedraw(base_size = 14)


###################################
# nMDS pigments
###################################

res.nmds <- metaMDS(PERCENT,"bray",k=2)
summary(res.nmds)

ggplot(data.frame(res.nmds$points),
       aes(x=MDS1,y=MDS2,col=fac))+
  geom_point(aes(size=lab))+
  stat_ellipse()+
  theme_linedraw(base_size = 14)





