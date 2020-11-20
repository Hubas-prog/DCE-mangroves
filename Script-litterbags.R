###################################
# DCE Mangroves
###################################

###################################
# Packages
###################################
library(vegan)
library(ade4)
library(FactoMineR)
library(wordcloud)
library(ggplot2)

#============================ ICI C'EST CONCARNEAU======

percent=read.table("/Users/cedric.hubas/Documents/CEDsync/DATASETS A VALORISER/GT DCE/Camille Pigments/pourcentages.txt",h=T)
PERCENT=read.table("/Users/cedric.hubas/Documents/CEDsync/DATASETS A VALORISER/GT DCE/Camille Pigments/POURCENT.txt",h=T)
FA=read.table("//Users/cedric.hubas/Documents/CEDsync/DATASETS A VALORISER/GT DCE/Baptiste FA/tables/FA.table.percent.txt",h=T)
FA$rep=substr(FA$groupe,3,4)
FA$station=substr(FA$groupe,1,2)
FA$time=substr(FA$groupe,5,6)
FA$ID=as.factor(paste(FA$rep,FA$station,FA$time,sep=""))

plot(percent[,1]~percent[,2])
names(percent)
pairs(percent)
plot(percent$Pheophytine.a ~percent$Pheophytine.a1)
res=lm(percent$Pheophytine.a ~percent$Pheophytine.a1)
summary(res)
abline(res)

ggplot(percent,aes(y=Pheophytine.a,x= Pheophytine.a1)) + geom_point() + geom_smooth(method="lm",aes(col="red",fill="red")) + theme_bw()


#=======ACP-pigments==============
names(PERCENT)
multi=dudi.pca(PERCENT,scannf=F,nf=3)
par(mfrow=c(2,2))
varexp=multi$eig/sum(multi$eig)*100
plot(multi$li[,1:2],type="n",xlab=paste("axis1: ",round(varexp[1],2),"%"),ylab=paste("axis2: ",round(varexp[2],2),"%"))
fac=as.factor(paste(percent$station,percent$time))
points(multi$li[,1],multi$li[,2],pch=21,bg=c("red","red1","red2","red3","green","green1","green2","green3","blue","blue1","blue2","blue3")[fac],col="white",cex=2)
s.class(multi$li,fac,add.plot=T,axesell=F,cellipse=0,cpoint=0,clab=0.7,col=c("red","red1","red2","red3","green","green1","green2","green3","blue","blue1","blue2","blue3"))
plot(multi$co[,1:2],type="n",xlab=paste("axis1: ",round(varexp[1],2),"%"),ylab=paste("axis2: ",round(varexp[2],2),"%"),xlim=c(-1,1),ylim=c(-1,1))
#s.corcircle(multi$co,grid=F)
arrows(0,0,multi$co[,1],multi$co[,2],length = 0.1, angle = 30,col="grey")
text(multi$co[,1],multi$co[,2],names(PERCENT),col="blue",cex=0.8)
abline(h=0,lty=2)
abline(v=0,lty=2)

plot(multi$li[,c(1,3)],type="n",xlab=paste("axis1: ",round(varexp[1],2),"%"),ylab=paste("axis3: ",round(varexp[3],2),"%"))
points(multi$li[,1],multi$li[,3],pch=21,bg=c("red","red1","red2","red3","green","green1","green2","green3","blue","blue1","blue2","blue3")[fac],col="white",cex=2)
s.class(multi$li[,c(1,3)],fac,add.plot=T,axesell=F,cellipse=0,cpoint=0,clab=0.7,col=c("red","red1","red2","red3","green","green1","green2","green3","blue","blue1","blue2","blue3"))
plot(multi$co[,c(1,3)],type="n",xlab=paste("axis1: ",round(varexp[1],2),"%"),ylab=paste("axis3: ",round(varexp[3],2),"%"),xlim=c(-1,1),ylim=c(-1,1))
#s.corcircle(multi$co,grid=F)
arrows(0,0,multi$co[,1],multi$co[,3],length = 0.1, angle = 30,col="grey")
text(multi$co[,1],multi$co[,3],names(PERCENT),col="blue",cex=0.8)
abline(h=0,lty=2)
abline(v=0,lty=2)

exp=inertia.dudi(multi,col=T)
mean(exp$col.abs[,1]/100)

multi2=bca(multi,fac,scannf=F,nf=2)
plot(multi2)
multi2$ls
par(mfrow=c(1,2))
plot(multi2$ls[,c(1,2)],type="n",pch=21,bg=c("red","red1","red2","red3","green","green1","green2","green3","blue","blue1","blue2","blue3")[fac],col="white",cex=2)
s.class(multi2$ls,fac,add.plot=T,axesell=F,cellipse=0,cpoint=0,clab=0.7,col=c("red","red1","red2","red3","green","green1","green2","green3","blue","blue1","blue2","blue3"))
#plot(multi2$co[,c(1,2)],type="n",xlim=c(-1,1),ylim=c(-1,1))
s.corcircle(multi2$co,grid=F)


par(mfrow=c(1,3),las=2)
boxplot(PERCENT$Pha~percent$station*percent$time,border=c("red","green3","blue"),ylab="Ph?ophytines (%)")
legend("topleft",levels(percent$station),text.col=c("red","green3","blue"))
boxplot(PERCENT$Cb/PERCENT$Ca~percent$station*percent$time,border=c("red","green3","blue"),ylab="Rapport Cb:Ca")
legend("topleft",levels(percent$station),text.col=c("red","green3","blue"))
boxplot(PERCENT$L/PERCENT$Ca~percent$station*percent$time,border=c("red","green3","blue"),ylab="Rapport L:Ca")
legend("topleft",levels(percent$station),text.col=c("red","green3","blue"))


MY.DATA=data.frame(Pheophytines=PERCENT$Pha,Temps=c(5,10,21,31)[as.factor(percent$time)],Stations=percent$station)
ggplot(MY.DATA,aes(y= Pheophytines,x=Temps,col=Stations)) + geom_point() + geom_smooth(method="lm") + ylab("Ph?ophytine (%)") + xlab("Temps (j)")


M1=lm(MY.DATA$Pheophytines~ MY.DATA$Temps*MY.DATA$Stations)
anova(M1)
coef(M1)

M2=lm(MY.DATA$Pheophytines[MY.DATA$Stations!="S2"]~ MY.DATA$Temps[MY.DATA$Stations!="S2"]*MY.DATA$Stations[MY.DATA$Stations!="S2"])
anova(M2)


MY.DATA.S13= MY.DATA[MY.DATA$Stations!="S2",]
ggplot(MY.DATA.S13,aes(y= Pheophytines,x=Temps,col=Stations)) + geom_point() + geom_smooth(method="lm") + ylab("Ph?ophytine (%)") + xlab("Temps (j)")




names(FA)
boxplot(FA$"X17.1w9"/FA$"X17.1w7"~FA$station*FA$time,border=c("red","green3","blue"),ylab="Rapport L:Ca")
legend("topleft",levels(percent$station),text.col=c("red","green3","blue"))

PERCENT$ID=percent$ID
temps=substr(PERCENT$ID,5,6)
temps=gsub("T1",5,temps)
temps=gsub("T2",10,temps)
temps=gsub("T3",20,temps)
temps=gsub("T4",30,temps)
temps=as.numeric(temps)

dac=data.frame(PERCENT,time=temps,station=as.factor(substr(PERCENT$ID,3,4)))
names(dac)
ggplot(dac,aes(y=Ca/L,x=time,col=station)) + geom_point() + geom_smooth(method="lm") + ylab("Ca:L ratio") + xlab("time (days)") + theme_bw()

mod1=lm(dac$Ca/dac$L~dac$time*dac$station);summary(mod1)
mod2=lm(dac$Ca/dac$L~dac$time+dac$station);summary(mod2)
mod3=lm(dac$Ca/dac$L~dac$time);summary(mod3)
anova(mod1,mod2)
anova(mod2,mod3)
anova(mod1,mod3)

#S1-S3
anova(lm(dac$Pha[dac$station!="S2"]~dac$time[dac$station!="S2"]+dac$station[dac$station!="S2"]),lm(dac$Pha[dac$station!="S2"]~dac$time[dac$station!="S2"]))
#S1-S2
anova(lm(dac$Pha[dac$station!="S3"]~dac$time[dac$station!="S3"]+dac$station[dac$station!="S3"]),lm(dac$Pha[dac$station!="S3"]~dac$time[dac$station!="S3"]))
#S2-S3
anova(lm(dac$Pha[dac$station!="S1"]~dac$time[dac$station!="S1"]+dac$station[dac$station!="S1"]),lm(dac$Pha[dac$station!="S1"]~dac$time[dac$station!="S1"]))

#======= PCoA pig ==============

names(PERCENT)
multi=dudi.pco(cailliez(vegdist(PERCENT,"bray")),scannf=F,nf=3)

varexp=multi$eig/sum(multi$eig)*100
plot(multi$li[,1:2],type="n",xlab=paste("axis1: ",round(varexp[1],2),"%"),ylab=paste("axis2: ",round(varexp[2],2),"%"))
fac=as.factor(paste(percent$station,percent$time))
points(multi$li[,1],multi$li[,2],pch=21,bg=c("red","red1","red2","red3","green","green1","green2","green3","blue","blue1","blue2","blue3")[fac],col="white",cex=2)
s.class(multi$li,fac,add.plot=T,axesell=F,cellipse=0,cpoint=0,clab=0.7,col=c("red","red1","red2","red3","green","green1","green2","green3","blue","blue1","blue2","blue3"))
plot(multi$co[,1:2],type="n",xlab=paste("axis1: ",round(varexp[1],2),"%"),ylab=paste("axis2: ",round(varexp[2],2),"%"),xlim=c(-1,1),ylim=c(-1,1))
s.corcircle(multi$co,grid=F)
arrows(0,0,multi$co[,1],multi$co[,2],length = 0.1, angle = 30,col="grey")
text(multi$co[,1],multi$co[,2],names(PERCENT),col="blue",cex=0.8)
abline(h=0,lty=2)
abline(v=0,lty=2)

#======= double PCoA FA/pig ==============
names(PERCENT)
multi=dudi.pco(cailliez(vegdist(PERCENT,"bray")),scannf=F,nf=3)



plot(coinertia(ant1, gen1, scann = FALSE))




#=======ACP-AG==============
names(FA)
multi=dudi.pca(FA[,2:46],scannf=F,nf=2)
fac=as.factor(paste(FA$station,FA$time))
multi2=bca(multi,fac,scannf=F,nf=2)
plot(multi2)
multi2$ls
par(mfrow=c(1,2))
plot(multi2$ls[,c(1,2)],type="n",pch=21,bg=c("red","red1","red2","red3","green","green1","green2","green3","blue","blue1","blue2","blue3")[fac],col="white",cex=2)
s.class(multi2$ls,fac,add.plot=T,axesell=F,cellipse=0,cpoint=0,clab=0.7,col=c("red","red1","red2","red3","green","green1","green2","green3","blue","blue1","blue2","blue3"))
#plot(multi2$co[,c(1,2)],type="n",xlim=c(-1,1),ylim=c(-1,1))
s.corcircle(multi2$co,grid=F)


#============================ ICI C'EST LA STATION : 
#========= MFA totale pigments/ FA ======

PERCENT$ID=percent$ID
FA$ID

DATA=merge(PERCENT,FA,by="ID")
names(DATA)
DATA2=DATA[,-c(1,14,60:62)]
names(DATA2)
colnames(DATA2)=gsub("X","",names(DATA2))
colnames(DATA2)=gsub(".",":",names(DATA2),fix=T)

write.table(DATA2)


mfa.res=MFA(DATA2, group=c(12,45), type=c("s","s"),name.group=c("pigments","acides gras"),graph=F)
plot(mfa.res$ind$coord[,1],mfa.res$ind$coord[,2],type="p",pch=21,bg=rainbow(3)[as.factor(DATA$station)],col="white",cex=1.5,xlab=paste(round(mfa.res$eig[1,2],2),"%"),ylab=paste(round(mfa.res$eig[2,2],2),"%"))
ordihull(mfa.res$ind$coord[,1:2],paste(DATA$time,DATA$station),lab=T,col="grey")

plot(mfa.res,choix="var",cex=0.7)



#======== MFA inter-groupes (pour mettre en ?vidence un facteur donn?)

DATA3=data.frame(DATA[,2:13]/dudi.pca(DATA[,2:13],scannf=F,nf=3)$eig[1],DATA[,15:59]/dudi.pca(DATA[,15:59],scannf=F,nf=3)$eig[1])

pca.res=dudi.pca(DATA3,scannf=F,nf=3)
names(DATA3)
biplot(pca.res)
inter=bca(pca.res,as.factor(DATA$station),scannf=F,nf=2)
plot(inter)
inter$ratio

par(mfrow=c(1,2))
plot(inter$ls,pch=21,cex=2,col="white",bg=c("purple","orange","green")[as.factor(DATA$station)])
ordihull(inter$ls,as.factor(paste(DATA$station,DATA$time)),lab=T)
plot(inter$co,xlim=c(-1,1),ylim=c(-1,1),type="n")
text=gsub("X","",rownames(inter$co))
text=gsub(".",":",text,fix=T)
arrows(0,0,inter$co[,1],inter$co[,2],col="grey",length=0.1)
textplot(inter$co[,1],inter$co[,2],text,col=c("red","blue")[as.factor(c(rep("pigments",12),rep("fatty acids",45)))],cex=0.7,new=F)
symbols(0,0,circle=1,add=T,inches=2)



