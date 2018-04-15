setwd("E:/Rwd2")
srcdir <- "./"
#### BOOT THOSE BAD BOIZ
library(bio3d)
library(plot3D)
#library(cairoDevice)

####Import Data for RMSD analysis
trjrmsd <- NULL
trjdat <- NULL
filedat <-NULL
for(i in c("5a1412")) { 
  for(j in 1:3500) {  
    filename <- paste(srcdir,"md_",i,"_",j,"-com.dat",sep="")
    print(filename)
    
    if(file.exists(filename)) {
      tmp <- read.table(filename,skip=2)
      #print(j)   
      filedat <- c(filedat,paste(i,"_",j,sep=""))
      # RMSD of ligand
      trjrmsd <- rbind(trjrmsd,tmp[,5])
      # X,Y,Z data on ligand CoM
      trjdat <- rbind(trjdat,unlist(tmp[,2:4]))
    }
  }
}
trjarr <- array(as.matrix(trjdat),dim=c(3500,251,3),dimnames=list(filedat,NULL,c("X","Y","Z")))
rownames(trjrmsd) <- filedat


#TRJDAT DATA RMSD ANALYSIS FOR SF
dim(trjrmsd)
str(trjrmsd)
head(trjrmsd)

#Convert the trjrmsd imported data file into a csv file that can be read in as a dataframe 
write.csv(trjrmsd, file="trjrmsd_tablesf.csv")
trjrmsd_table<-read.csv("trjrmsd_tablesf.csv")
trjrmsd_table$X<- c(1:3500)
trjrmsd2<-trjrmsd_table[,-1]

##Load in relevent data from Network Analysis
load("Clustercom Rigid Vac.RData")
load("clushub rigid vac.RData")
nclus=8
clushub <- (clushub*2) - 35
load("clustnumb Rigid Vac.RData")
#reprocess clustnumb just post the hard one, these are quick enough im not bothered 
#remove trajectories that dont go near any cluster points - this can be done in previous loop
names(clustnumb) <- seq_along(clustnumb)
clustnumb[sapply(clustnumb, is.null)] <- NULL

#assigns the trajectory points close to a cluster point to that clusters number
clustercom[,4] <- as.numeric(clustercom[,4])  
for (i in 1:length(clustnumb)) {
  clustpoint2 <- NULL
  clustpoint <- clustercom[rownames(clustercom) %in% clustnumb[[i]][,1],]
  for (y in 1:nrow(clustnumb[[i]])){
    c <- clustnumb[[i]][y,]
    c2 <- which(rownames(clustpoint) == c)
    clustpoint2 <- (rbind(clustpoint2, clustpoint$`kcluster$cluster`[c2]))
  }
  clustnumb[[i]] <- cbind(clustnumb[[i]], clustpoint2)
  percen <- as.integer((i/length(clustnumb)*100))
  print(percen)
}

#make cluster numer a factor
for (p in 1:length(clustnumb)) {
  clustnumb[[p]][,2] <- factor(clustnumb[[p]][,2])
  print(p)
}


###LOOKS FOR MULT_CLUSTERS 
mult_cluster <- NULL
for (b in 1:length(clustnumb)){
  factornumb <- table(clustnumb[[b]][,2])
  if (nrow(factornumb) >= 2) { #### a=nclus 
    mult_cluster <- c(mult_cluster, b)
    mult_cluster2 <- NULL
    clustnumbnames <- names(clustnumb)
  }
}
for (c in 1:length(mult_cluster)) {
  mult_cluster2 <- c(mult_cluster2, clustnumbnames[mult_cluster[c]])
  tmpsub <- clustnumb[mult_cluster2]
}

###SUBSET DATA BY IF THEY VISIT A CLUSTER
clus1trj <- NULL
clus2trj <- NULL
clus3trj <- NULL
clus4trj <- NULL
clus5trj <- NULL
clus6trj <- NULL
clus7trj <- NULL
clus8trj <- NULL
for (i in 1:length(clustnumb)){
  clu1tab <- as.data.frame(table(clustnumb[[i]][,2]))
  clu1boolean <- clu1tab[,1] == 1
  clu2boolean <- clu1tab[,1] == 2
  clu3boolean <- clu1tab[,1] == 3
  clu4boolean <- clu1tab[,1] == 4
  clu5boolean <- clu1tab[,1] == 5
  clu6boolean <- clu1tab[,1] == 6
  clu7boolean <- clu1tab[,1] == 7
  clu8boolean <- clu1tab[,1] == 8
  if (any(clu1boolean)) {
    
    clus1trj <- c(clus1trj, as.integer(names(clustnumb)[[i]]))
  }
  if(any(clu2boolean)) {
    
    clus2trj <- c(clus2trj, as.integer(names(clustnumb)[[i]]))
  }
  if(any(clu3boolean)) {
    
    clus3trj <- c(clus3trj, as.integer(names(clustnumb)[[i]]))
  }
  if(any(clu4boolean)) {
    
    clus4trj <- c(clus4trj, as.integer(names(clustnumb)[[i]]))
  }
  if(any(clu5boolean)) {
    
    clus5trj <- c(clus5trj, as.integer(names(clustnumb)[[i]]))
  }
  if(any(clu6boolean)) {
    
    clus6trj <- c(clus6trj, as.integer(names(clustnumb)[[i]]))
  }            
  if(any(clu7boolean)) {
    
    clus7trj <- c(clus7trj, as.integer(names(clustnumb)[[i]]))
  }
  if(any(clu8boolean)) {
    
    clus8trj <- c(clus8trj, as.integer(names(clustnumb)[[i]]))
  }
} 



###PROCESS DATA -  THE RMSD DATA 
meltout<-NULL
trjdiff <- NULL
hisout<-NULL
nclus<-8
colhis <- get_colors(c(1:nclus), group.col = brewer.pal(nclus, "Set2"))
for (i in 1:nclus){
  loopnam <- paste("clus",i,"trj",sep="")
  lp2 <- paste("trjrmsd_Retained",i,sep="")
  lp3 <- paste("avRMSDC",i,sep="")
  assign(lp2,trjrmsd2[get(loopnam),,drop=F])
  assign(lp3, cbind(1:251, as.data.frame(colMeans(get(lp2)))))
  FinalRMSD<-get(lp2)[,251]
  lpbind <- cbind(melt(FinalRMSD),i)
  meltout <- rbind(meltout,lpbind)
  trjlast<- get(lp2)[,ncol(get(lp2))]
  trjmin <- as.vector(apply(get(lp2),1,min))
  trjdiff<- trjmin-trjlast
  trjdiff<- cbind(melt(trjdiff),i)
  hisout<- rbind(hisout,trjdiff)
}

##PREP COLOUR
colhis <- hisout$i
colhis <- get_colors(colhis, group.col = brewer.pal(nclus, "Set2"))
collis <- get_colors(c(1:8),group.col = brewer.pal(nclus, "Set2") )
##STACKED HISTOGRAM of min value of the trajectory vs final - not used in Report as useless

colnames(meltout) <- c("RMSD","Basin")
#PREP COLOURS
colfin <- meltout$Basin
colfin <- get_colors(colfin, group.col = brewer.pal(nclus, "Set2"))


###AVERAGED RMSD - rename ####
colnames(avRMSDC1) <- c("Basin 1"," ")
colnames(avRMSDC2) <- c("Basin 2"," ")
colnames(avRMSDC3) <- c("Basin 3"," ")
colnames(avRMSDC4) <- c("Basin 4"," ")
colnames(avRMSDC5) <- c("Basin 5"," ")
colnames(avRMSDC6) <- c("Basin 6"," ")
colnames(avRMSDC7) <- c("Basin 7"," ")
colnames(avRMSDC8) <- c("Basin 8"," ")


RMSDplot1<- ggplot(avRMSDC1, aes(`Basin 1`,` `)) + geom_line()+
  theme_classic()+
  ylim(c(15, 50))
RMSDplot2<- ggplot(avRMSDC2, aes(`Basin 2`,` `)) + geom_line()+
  theme_classic()+
  ylim(c(15, 50))
RMSDplot3<- ggplot(avRMSDC3, aes(`Basin 3`,` `)) + geom_line()+
  theme_classic()+
  ylim(c(15, 50))
RMSDplot4<- ggplot(avRMSDC4, aes(`Basin 4`,` `)) + geom_line()+
  theme_classic()+
  ylim(c(15, 50))
RMSDplot5<- ggplot(avRMSDC5, aes(`Basin 5`,` `)) + geom_line()+
  theme_classic()+
  ylim(c(15, 50))
RMSDplot6<- ggplot(avRMSDC6, aes(`Basin 6`,` `)) + geom_line()+
  theme_classic()+
  ylim(c(15, 50))
RMSDplot7<- ggplot(avRMSDC7, aes(`Basin 7`,` `)) + geom_line()+
  theme_classic()+
  ylim(c(15, 50))
RMSDplot8<- ggplot(avRMSDC8, aes(`Basin 8`,` `)) + geom_line()+
  theme_classic()+
  ylim(c(15, 50))

###Plot Appendix 7.5
grid.draw(rbind(ggplotGrob(RMSDplot1), ggplotGrob(RMSDplot3), ggplotGrob(RMSDplot5), ggplotGrob(RMSDplot2), size = "last"))
grid.draw(rbind(ggplotGrob(RMSDplot6), ggplotGrob(RMSDplot7), ggplotGrob(RMSDplot8), ggplotGrob(RMSDplot4), size = "last"))

#####Plot Appendix 7.4 - TYPICAL RMSD - 1 mobile as 3 clus 
typ <- trjrmsd2[1,]
typ <- as.data.frame(t(typ))
dim(typ)
typ <- cbind(c(1:251), typ)
colnames(typ) <- c("Time","DHF RMSD away from DHFR (Å)")
View(typ)
typplot<- ggplot(typ, aes(Time,`DHF RMSD away from DHFR (Å)`)) + geom_line()+
  theme_classic()+
  ylim(c(15, 50))+
  ylab("DHF RMSD away from DHFR (Å)")
typplot
