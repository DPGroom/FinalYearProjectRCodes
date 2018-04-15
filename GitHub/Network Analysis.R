library('RColorBrewer')

#FOR CLUSTERING
library(multcomp)
library(kml3d)

## FOR PLOTTING
library(rgl)
library(bio3d)
library(ggfortify)
library(ggplot2)

###Network
library(igraph) 
library(markovchain)
require(reshape2)

#
library(misc3d)
library(sm)
library(fields)
library(rdist) ## defo this one





#### FUNCTIONS ####
###get a decent legend that can be added to a grid layout, taken from https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 




##FORCE DP taken from https://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

#get rep row in taken from https://www.r-bloggers.com/a-quick-way-to-do-row-repeat-and-col-repeat-rep-row-rep-col/
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

###taken from http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}




####LOAD IN TO GET DATA IN CORRECTLY FOR PROCESSING #####
setwd("G:/")

####First import the data
srcdir <- "./"
# md_5a1412_1-com.dat
trjdatxyz <- NULL
runnumb <- 1
for(i in c("5a1412")) { 
  for(j in 1:3500) {  
    filename <- paste(srcdir,"md_",i,"_",j,"-com.dat",sep="")
    print(filename)
    
    if(file.exists(filename)) {
      #print("debug")
      tmp <- read.table(filename,skip=2, col.names = c("Number", "X", "Y", "Z", "RMSD"))
      tmp2 <- tmp[,2:4]
      tmp3 <- cbind(tmp2, runnumb)
      runnumb <- runnumb + 1
      trjdatxyz <- rbind(trjdatxyz,tmp3)
    }
    
  }
  
}




####PROTEIN - Plotting/Import####
#IMPORT AND SELECT CA
#set the range for plotting (rgl is inconsistent at sticking to these limits)
rg <- c(-60,60)
#first read in the initial crd, only the protein is considered, the ligand is ignored
crd <- read.crd.charmm("md_5a1412_1_init.crd",ext=T)
ca.selnope <- which(crd$calpha)
##an atom within the Ligand was also selected so removed
ca.sel <- ca.selnope[1:186]
xyz <- matrix(crd$xyz[atom2xyz(ca.sel)],ncol=3,byrow=T)



#Converts from data.frame to array 
trjarr2 <- simplify2array(by(trjdatxyz[1:3], trjdatxyz$runnumb, as.matrix))
save(trjarr2, file="trajarr2.RData")

#SET LIMITS for the KDE - Methods 2.2.1
limnu<- 35
slim<-c(-limnu,limnu)
boxlim<-c(slim,slim,slim)


#DO THE ESTIMATE 
attach(trjdatxyz)
Trj3D<-kde3d(X,Y,Z, n=limnu, lims=boxlim)
#save(Trj3D,file="Trj3D.RData")
#load("Trj3D Rigid Vac.RData")


#multiply kernel density data so it is easier to work with
lev<-Trj3D$d*1000



#shows cluster number above a chosen level. (this probably need changing depending on your data)
hist(lev[lev>0.1])
levv <- 0.15

##Plot the KDE  - Figure 3.3A
lines3d(xyz[,1],xyz[,2],xyz[,3],col="black", xlim=rg,ylim=rg,zlim=rg,phi =30,theta =30, box=T, lit=F)+
  contour3d(lev,level=levv,Trj3D$x,Trj3D$y,Trj3D$z, color= "#08306b", color2= "white", add=T, engine="rgl", lit=T)+
  axes3d(edges = "bbox", labels = TRUE, tick = TRUE, nticks = 5, box = FALSE, expand = 1.5, xlab="x",ylab="y",zlab="z")

#export
rgl.postscript("../", fmt = "pdf", drawText = TRUE )


###STARTS CLUSTERING BY KDE to resolve Basins ####
#finds kernel density cells above a certain level. assigns TRUE or FALSE
d.data <- lev
signif <- d.data >= levv
nclus<- 8

#Finds significant cluster points. Co-ords calculated from these points. (depending on n and lims value set earlier)
cluster <- which(signif, arr.ind = TRUE)
cluster_df <- as.data.frame(cluster)
clustercom <- as.data.frame(cluster*2- limnu)

#Run k means clustering on dist data
kcluster <- kmeans(cluster_df, nclus, nstart = 20)
kcluster$cluster <- as.factor(kcluster$cluster)

#add cluster number calculated from kmeans to cluster point
clustercom <- cbind(clustercom, kcluster$cluster)
#save(clustercom, file="Clustercom.RData")
#load("Clustercom Rigid Vac.RData")

kmeansc <- clustercom$`kcluster$cluster`
clushub <- kcluster$centers
#save(clushub, file="clushub.RData")
#load("clushub rigid vac.RData")

dimcom <- dim(clushub)
dimcom <- dimcom[1]
## HUBS = CLUSTER CENTRES 
## CLUSOINTS = CLUSTER POINTS 
#translate
clushub <- (clushub*2) - 35
dimhub <- dim(clushub)
dimhub <- dimhub[1]  


####Principal Coordinate Analysis of cluster points####### - Methods 2.2.1 - Groom 2018
#make a copy just in case
clustercom2 <- cbind(clustercom)
colnames(clustercom2) <- c("dim1","dim2","dim3","Cluster")
###APPEND CENTRES for accurate labels
Clushubb <- cbind(clushub,Cluster=1:nclus)
Cluscom_w_hub <- rbind(Clushubb, clustercom2)
#pairwise distances
dis_cluscom_w_HUB<- cdist(Cluscom_w_hub[1:3],Cluscom_w_hub[1:3])
#the analysis
PC_CwH1 <- prcomp(dis_cluscom_w_HUB,center = TRUE,scale. = TRUE)
#save(PCi_CwH1, file="PCi_wH1 Rigid Vac.RData")
#load("PCi_wH Rigid Vac.RData")



#extract x/y for ggplot
PCi_CwH <- data.frame(PC_CwH1$x, Cluster=Cluscom_w_hub$Cluster)
#save(PCi_CwH, file="PCi_wH Rigid Vac.RData")
#or alt. use autoplot, limited aes commands
autoplot(PC_CwH1)


###GGPLOT PCoA
###Plot Topographic  - Figure 3.3C
dtmp <- dim(PCi_CwH)[1]
plot2 <- ggplot(PCi_CwH[(nclus+1):dtmp,],aes(x=PC1,y=PC2, col=factor(Cluster)), show.legend=FALSE)+
  geom_point(alpha=0.5, show.legend=FALSE)+
  theme_classic()+
  labs(colour="Cluster")+
  scale_colour_brewer(palette = "Set2", direction = 1)+
  theme(legend.box.background = element_rect(),legend.box.margin = margin(4, 4, 4, 4))+
  geom_point(data=PCi_CwH[1:nclus,],aes(x=PC1,y=PC2, col=factor(Cluster)), size=7,alpha= 0.9, show.legend=T)+
  geom_text(data=PCi_CwH[1:nclus,],aes(label=Cluster, fontface= "bold"), col= "grey 3", show.legend = FALSE)+
  expand_limits(x=c(-10,15), y=c(-8,8)) + xlab("PC1 = 76.54%") + ylab("PC2 = 16.10%")+ labs(colour="Basin")
plot2
legend <- g_legend(plot2) 
grid.draw(legend) 

##FOR an faded circle around centres, not used within study
geom_point(data=PCi_CwH[1:nclus,],aes(x=PC1,y=PC2, col=factor(Cluster)), size=8,alpha= 0.2, show.legend=FALSE)+


  
###Assigning Trajectories to Basins - - Methods 2.2.2 - Thanks to Tinsley (2018)
detach()
attach(clustercom)
start.time <- Sys.time()
clustnumb <- list()
trjdim <- dim(trjarr2)
trjnumber <- trjdim[3]
for (i in 1:trjnumber){
  trjrun <- trjarr2[,,i]
  disrun <- cdist(trjrun, clustercom[,1:3])
  coldisrun <- as.character(seq(1:ncol(disrun)))
  dfdisrun <- as.data.frame(disrun)
  colnames(dfdisrun) <- coldisrun
  dfcolname <- NULL
  for (y in 1:nrow(dfdisrun)){
    if (min(dfdisrun[y,]) <= 2) {
      colname<- as.data.frame(names(dfdisrun)[which.min(apply(dfdisrun[y,],MARGIN=2,min))], stringsAsFactors = FALSE)
      rownames(colname) <- c(y)
      colnames(colname) <- c("Cluster Column Number")
      dfcolname <- rbind(dfcolname, colname) 
    }
    clustnumb[[i]] <- dfcolname
  }
  percen <- as.integer((i/trjnumber*100))
  print(percen)
}
end.time <- Sys.time()
end.time-start.time 
#save(clustnumb, file="clustnumb Rigid Vac.RData")
#load("clustnumb Rigid Vac.RData")


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


###LOOKS FOR MULT_CLUSTERS - Methods 2.2.2
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



###CONSTRUCT NETWORK - GOLDEN EGG PATCH - Methods 2.2.3 - Groom (2018)

# This loop only samples mobile trajectories unfortunately
# This is the final loop, earlier iterations have been excluded
tempnet <- NULL
output<- NULL
##SAMPLE MULTICLUSTERS FOR RELATIONAL INFO/EDGE/PROP - Methods 2.2.3
for (d in clustnumbnames) {
  tmpsubb <- clustnumb[[d]][["clustpoint2"]]
  rletmp <- rle(as.numeric(levels(tmpsubb)[as.integer(tmpsubb)]))
  ###This is a list of the basins for each traj.
  hubnet<- rletmp$values
  ###This is the length of runs 
  hubkin <- rletmp$lengths
  ##Looks at how long the run is, uses that to make sure a traj. is all accounted for
  ln <- dim(as.data.frame(hubnet))[1]
  tempd2 <- NULL
  d2pr <-NULL
  for (e in 1:ln) {
    ###check to see if e = penultimate migration, if it's the final one else loop won't work/nomigration in final value thus no edge
    if ((e+1 > ln)==TRUE) {
      warningg <- paste(d,"has been processed, calculating probability", sep=" ")
      print(warningg)
    }
    #takes the basin in the rle output, identifies adjacent basin, records the length of time spent within the original basin
    else { 
      tempnet <-data.frame(hubnet[e], hubnet[e+1],hubkin[e])
      colnames(tempnet) <- c("N1","N2","Time")
      d2pr <- rbind(d2pr, tempnet)
    }
  }
  for (z in 1:nclus){
    if (any(d2pr$N1 == z)==FALSE){
    }
    ###Early adaptations involved simply appending on each trajectory to a master file, rbind ran into issues with RAM; aggregation serves as a compression method
    else {
      subset1<- d2pr[which(d2pr$N1==z),,]
      aggr <- aggregate(subset1, by=list(subset1$N2), FUN=sum)
      aggr2 <- as.data.frame(cbind(z,aggr$Group.1,(aggr$N1/z),aggr$Time))
      colnames(aggr2) <- c("N1","N2","nmigrate","Time")
      output <- rbind(output, as.data.frame(aggr2))
      colnames(output) <- c("N1","N2","nmigrate","Time")
    }
  }
}
###Initial Output is all of the aggr2's which is then reaggregated to give an accurate prop of the whole population



###reaggregates the output in a simialr method to Orig but results in a final graph
W_edgelist <- NULL
for (i in 1:nclus){
  subset2<- output[which(output$N1==i),,]
  aggrout <- aggregate(subset2, by=list(subset2$N2), FUN=sum)
  ##calc prob of retention
  mrkout <- 1-(sum(aggrout$nmigrate)/sum(aggrout$Time))
  mrkpout <- as.vector(c(i,i,mrkout))
  #calc prob of leaving
  pbout <- aggrout$nmigrate/sum(aggrout$Time)
  elout <- cbind(i,aggrout$Group.1,pbout)
  elout <- rbind(elout,as.numeric(mrkpout))
  W_edgelist<- rbind(W_edgelist,elout)
}




## Give them nice names 
colnames(W_edgelist) <- c("N1","N2","Probability")
write.csv(W_edgelist, file="Weighted Edge List - Rigid Vac.csv")
#W_edgelist <- read.csv("Weighted Edge List - Rigid Vac.csv",sep=",")
#W_edgelist$X<-NULL


##For Figure 3.4 - simplify values with a prob of >=0.005
wammend2 <- as.data.frame(W_edgelist[W_edgelist[,3]>=0.005,,drop=F])
wammend2$Probability <- specify_decimal(wammend2$Probability,3)


###for plotting later on sep the markov from not
wammendmark <- wammend2[wammend2$N1==wammend2$N2,,drop=F]
wammendsubset <- wammend2[wammend2$N1!=wammend2$N2,,drop=F]




###PLOT THE iGRAPH
net<- graph_from_data_frame(wammend2)

#PREP COLOURS
colr <- wammend2$N1
colr<- get_colors(colr, group.col = brewer.pal(nclus, "Set2"))

####NON Topological - Methods 2.2.5 - Figure 3.4A/B
plot(net, edge.label=wammend2$Probability, 
     vertex.color=brewer.pal(nclus, "Set2"),vertex.label.color="Black", 
     edge.curved=TRUE, edge.arrow.size=.3, edge.label.color="black", edge.color=colr, asp=0.8)


###TOPGRAPHIC representation of network - for DHFR this isn't as good but functionality is there
nodes.coord <-PCi_CwH[1:nclus,1:2]
plot(net, rescale=FALSE, axes=TRUE, layout=as.matrix(nodes.coord[,c("PC1","PC2")]), vertex.size= 55,
     edge.label=wammend2$Probability, vertex.color=brewer.pal(nclus, "RdYlBu"),
     vertex.label.color="black", edge.curved=seq(-1.1, 1.1, length = ecount(net)),edge.arrow.size=.4,edge.label.color="black", edge.color=colr,
     edge.loop.angle=9.8, xlab="PC1 = 82.24%", ylab="PC2 = 10.29%",edge.label.x=-4.2,
     xlim=c(min(nodes.coord[,"PC1"])-2,max(nodes.coord[,"PC1"])+2), ylim=c(min(nodes.coord[,"PC2"])-2,max(nodes.coord[,"PC2"])+2),
     asp=0)




#### Plotting CLusterpoints:
###PREP COLOURS
colrr <- kmeansc
colrr<- get_colors(colrr, group.col = brewer.pal(nclus, "Set2"))


#Figure 3.3B
lines3d(xyz[,1],xyz[,2],xyz[,3],col="black", xlim=rg,ylim=rg,zlim=rg,phi =30,theta =30, box=T)+
  points3d(clustercom[,1], clustercom[,2], clustercom[,3], color = colrr,size=5)+
  axes3d(edges = "bbox", labels = TRUE, tick = TRUE, nticks = 5, box = F, expand = 1.5, xlab="x",ylab="y",zlab="z")
  view3d(userMatrix = newv)
rgl.postscript("../", fmt = "pdf", drawText = TRUE )

######### Intialising Colours
colrrsp <- c(1:nclus)
colrrsp <- get_colors(colrrsp, group.col = brewer.pal(nclus, "Set2"))

###initial plot for clusterpoints/centres##
intial <- function() {
  points3d(clustercom[,1], clustercom[,2], clustercom[,3], color = colrr) +
  rgl.spheres(clushub[,1],clushub[,2],clushub[,3], color = colrrsp, lit=F)
} 


########## ARROWS
###PREP COLOURS
colar <- wammendsubset$N1 
colar <- get_colors(colar, group.col = brewer.pal(nclus, "Set2"))
colnuc <- c(1:nclus)
colmark <- get_colors(colnuc, group.col = brewer.pal(nclus, "Set2"))


##PLOT THE ARROWS
arrow2 <- function(){
  for (i in 1:dim(wammendsubset)[1]){
    arrow3d(p0=clushub[wammendsubset$N1[i],], p1=clushub[wammendsubset$N2[i],],type="rotation", 
            n= 20, width=1/7, thickness=1/7, color=as.matrix(colar[i]), lit=FALSE, alpha= .5)
  }
}

##ACTUAL PLOTTING on a protein backbone  - Figure 4.1C
lines3d(xyz[,1],xyz[,2],xyz[,3],col="black", xlim=rg,ylim=rg,zlim=rg,phi =30,theta =30, box=T)+ 
axes3d(edges = "bbox", labels = TRUE, tick = TRUE, nticks = 5, box = F, expand = 1.5, xlab="x",ylab="y",zlab="z")+ 
intial()+
arrow2()+
points3d(clushub[wammendmark$N1,1],clushub[wammendmark$N1,2],clushub[wammendmark$N1,3], 
         color = as.matrix(colmark[wammendmark$N1]), lit=F, size=30, alpha=.5)
rgl.postscript("../", fmt = "pdf", drawText = TRUE )





#######Markov analysis##### - Methods 2.2.4
#load("Weighted Edge List - Rigid Vac.RData")
net2 <- graph_from_data_frame(W_edgelist)
ma_transitionProbs <- prop.table(as.matrix(get.adjacency(net2, attr="Probability")),1)
W_adjacency<-as.matrix(get.adjacency(net2, attr="Probability")) 
#write.csv(W_adjacency, file="Transition Matrix.csv")


tmA <- ma_transitionProbs
netmark <- new("markovchain",transitionMatrix= tmA ,
               name="RigidVacMarkov")
write.csv(netmark[,][,], file="Transition Matrix MARK.csv")
netmm <- netmark[,][,]
#save(netmm, file="Transition Matrix MARK.RData")

##Identify subgrapghs
communicatingClasses(netmark)



####STEADY STATE GRAPHS - Figure 3.6
arrayLop <- NULL
arwahay<-NULL
###genn= generation number, the max raised 'n' 
genn <- 50

##makes an array where each layer is the powered transition matrix
for (i in 1:genn){
  loopnam <- paste("T",i,sep="")
  nmloop <- netmark^i
  arrayLop <- cbind(arrayLop, as.matrix(nmloop[,][,]))
}
###NEED AN A WAY TO SUBSET THE ARRAY AND PLOT THE POINTS OVER TIME
arwahay <- array(data=arrayLop, dim = c( nclus , nclus , genn ))
arwahay[,,][,,]

##Generate all the plots in a loop
get_colors(c(1:nclus), group.col = brewer.pal(nclus, "Set2"))


df<-NULL
###DO A RUN WITH LEGEND to save the legend then redo with show.legend =F - plots each steady state graph
for (i in 1:nclus){
  dat<- as.data.frame(t(arwahay[,,][,i,1:genn]))
  colnames(dat) <- c(1:nclus)
  time<- c(1:genn)
  dat2 <- cbind(time, dat)
  df <- melt(dat2 ,  id.vars = 'time', variable.name = 'Cluster')
  finp <- df[df[1]==genn,,drop=F]
  lab <- cbind(genn,mean(finp$value[finp$value!=0]))
  df[df == 0] <- NA
  namaa<- paste("Cluster",i,sep = " ")
  lplot <- ggplot(df, aes(time,value)) +
    geom_line(aes(colour = Cluster),size=0.9, show.legend = F)+
    geom_point(aes(colour = Cluster),size=3,show.legend =F)+
    theme_classic()+
    labs(colour="Cluster")+
    scale_colour_brewer(palette = "Set2", direction = 1)+
    ylim(0,1) + xlab(namaa)+ ylab("")+  scale_x_continuous(breaks = c(5,10,15,20,25))+
    theme(legend.box.background = element_rect(),legend.box.margin = margin(4, 4, 4, 4))+
    annotate("text",label=specify_decimal(lab[2],3), x=24, y=(lab[2]+0.25))
  lpnam <- paste("TFPlot",i,sep="")
  assign(lpnam, lplot)
}

##save legend
legend <- g_legend(TFPlot1) 
grid.draw(legend) 



####THIS BIT IS BESPOKE FOR COMS that start at 0, ie if you were to run TFplot2 now it would miss n=1. Allows to seperate the communities
dat<- as.data.frame(t(arwahay[,,][,2,1:genn]))
colnames(dat) <- c(1:nclus)
time<- c(1:genn)
dat2 <- cbind(time, dat)
df <- melt(dat2 ,  id.vars = 'time', variable.name = 'Cluster')
df2 <- df[26:50,]
df4 <- df[76:100,]
df <- rbind(df2,df4)
colrr <- get_colors(c(1:nclus), group.col = brewer.pal(nclus, "Set2"))
finp <- df[df[1]==genn,,drop=F]
lab <- cbind(genn,mean(finp$value[finp$value!=0]))
TFPlot2 <- ggplot(df, aes(time,value)) +
  geom_line(aes(colour = df$Cluster),size=0.9, show.legend = F)+
  geom_point(aes(colour = Cluster),size=3,show.legend = F)+
  theme_classic()+
  labs(colour="Cluster")+
  scale_colour_manual(values=colrr[c("2","4")])+
  ylim(0,1) + xlab("Cluster 2")+ ylab("")+  scale_x_continuous(breaks = c(5,10,15,20,25))+
  annotate("text",label=specify_decimal(lab[2],3), x=24, y=(lab[2]+0.25))



###NICE Grid plots:) - Figure 3.5
grid.draw(rbind(ggplotGrob(TFPlot1), ggplotGrob(TFPlot3), ggplotGrob(TFPlot5), ggplotGrob(TFPlot2), size = "last"))
grid.draw(rbind(ggplotGrob(TFPlot6), ggplotGrob(TFPlot7), ggplotGrob(TFPlot8), ggplotGrob(TFPlot4), size = "last"))

       