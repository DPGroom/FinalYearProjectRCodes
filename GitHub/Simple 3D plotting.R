
#FOR CLUSTERING
library(multcomp)
library(kml3d)

#FOR PLOTTING PROTEIN ECT 

#LOWER CASE d 
library(bio3d)
library(plot3D)
library(rgl)

#CENTrE OF MASS
library(misc3d)
library(sm)
library(fields)
library(gplots)

##Network



#LOAD STUFF IN
setwd("G:/Rwd2")
####SF TEST DATA IN USB
srcdir <- "./"
# md_5a1412_1-com.dat
trjrmsd <- NULL
trjdat <- NULL
trjdatxyz <-NULL 
filedat <-NULL
filedat4xyz<- NULL
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
      trjdatxyz<-rbind(trjdatxyz,tmp[,2:4])
    }
  }
}
trjarr <- array(as.matrix(trjdat),dim=c(3500,251,3),dimnames=list(filedat,NULL,c("X","Y","Z")))
rownames(trjrmsd) <- filedat
rownames(trjdatxyz) <- filedat4xyz


####### ASSIGNING DATA FRAMESA ######
#FINAL FRAME
endpt <- dim(trjrmsd)[2]

#DEFINE CUT OFF
rmsdlim <- 30
gone <- unique(which(trjrmsd[,endpt] > rmsdlim,arr.in=T))

# reduce the data matrix, by removing the escapees
td <- trjarr[-gone,,,drop=F]
tdid <- rownames(trjarr)[-gone]
dim(td)[1]
#FINAL COORDINATE POSITION
tds <- matrix(td[1:dim(td)[1],endpt,1:3],nrow=dim(td)[1], ncol=3,byrow=F)

#Initial
tdin <- matrix(td[1:dim(td)[1],1,1:3],nrow=dim(td)[1],ncol=3,byrow=F)

#td as a matrix
tdmat <- matrix(td[1:dim(td)[1],1:251,1:3],nrow=dim(td)[1],ncol=3,byrow=F)


##### CLUSTERING#######
#Generate a correlation-based distance matrix OF FINAL COORDINATES 
dm<-1-cor(t(tds))
nclus<-8

#Generate K means cluster graphs with 4 (nclus) centers: 
km<-kmeans(tds,centers=nclus)

##ALL##
allidclu <- cbind(idclu$Cluster,tdmat)

##FINAL####
idclu <- cbind(km$cluster,tds)

##INITIAL###
inidclu<- cbind(km$cluster,tdin)




#WRITES OUT THE CLUSTER OF THE SUBSET
colnames(idclu) <- c("Cluster","X", "Y","Z")
write.csv(idclu, file="5a1214_Clusterd.csv")
write.csv(inidclu, file="5a1214_ClusterdINITIL.csv")
#idclu<-read.csv("5a1214_Clusterd.csv")
#idclu$X.1<-NULL
#inidclu <- read.csv("5a1214_ClusterdINITIL.csv")
#inidclu$X.1<-NULL


#LABELS THE CLUSTER
colnames(km$centers) <- c("X", "Y","Z")
write.csv(km$centers, file="5a1215_sfClusterSPOT.csv")
#ctre <- read.csv("5a1215_sfClusterSPOT.csv")
#ctre$X.1<-NULL
#ctre <- km$centers


#### GROUPING BY CLUSTER ##### THIS CAN BE IMPROVED MASSIVELY WITH LOOPING 



####ALL####
allidcludat <- as.data.frame((allidclu))
clutrj1 <- td[allidcludat[1]=="1",,,drop=F]
clutrj2 <- td[allidcludat[1]=="2",,,drop=F]
clutrj3 <- td[allidcludat[1]=="3",,,drop=F]
clutrj4 <- td[allidcludat[1]=="4",,,drop=F]
clutrj5 <- td[allidcludat[1]=="5",,,drop=F]
clutrj6 <- td[allidcludat[1]=="6",,,drop=F]
clutrj7 <- td[allidcludat[1]=="7",,,drop=F]
clutrj8 <- td[allidcludat[1]=="8",,,drop=F]


######FINAL####
#TO BE ABLE TO PROCESS INTO SUBSETS
idcludat <- as.data.frame((idclu))

#SUBSETTED, LOOP AT THE BOTTOM OF CODE NEEDS WORK
cluf1 <- idcludat[which(idcludat$Cluster=="1"),,]
cluf2 <- idcludat[which(idcludat$Cluster=="2"),,]
cluf3 <- idcludat[which(idcludat$Cluster=="3"),,]
cluf4 <- idcludat[which(idcludat$Cluster=="4"),,]
cluf5 <- idcludat[which(idcludat$Cluster=="5"),,]
cluf6 <- idcludat[which(idcludat$Cluster=="6"),,]
cluf7 <- idcludat[which(idcludat$Cluster=="7"),,]
cluf8 <- idcludat[which(idcludat$Cluster=="8"),,]


###INTIAL#####

##PROCESSING FOR SUBSETTING###
inidcludat <- as.data.frame((inidclu))
inidcludat$X<-NULL

#SUBSET, LOOP AT BOTTOW NEEDS SOME WORK

cluin1 <- inidcludat[which(idcludat$Cluster=="1"),,]
cluin2 <- inidcludat[which(idcludat$Cluster=="2"),,]
cluin3 <- inidcludat[which(idcludat$Cluster=="3"),,]
cluin4 <- inidcludat[which(idcludat$Cluster=="4"),,]
cluin5 <- inidcludat[which(idcludat$Cluster=="5"),,]
cluin6 <- inidcludat[which(idcludat$Cluster=="6"),,]
cluin7 <- inidcludat[which(idcludat$Cluster=="7"),,]
cluin8 <- inidcludat[which(idcludat$Cluster=="8"),,]


####################PLOTTING########################

#SET RANGE OF FIGURE
rg <- range(td)
rg <- c(-60,60)

####THE PROTEIN####

#IMPORT AND SELECT CA
crd <- read.crd.charmm("md_5a1412_1_init.crd",ext=T)

ca.selnope <- which(crd$calpha)
ca.sel <- ca.selnope[1:186]

xyz <- matrix(crd$xyz[atom2xyz(ca.sel)],ncol=3,byrow=T)






# plot the protein structure
lines3d(xyz[,1],xyz[,2],xyz[,3],col="black", xlim=rg,ylim=rg,zlim=rg,phi =30,theta =30, box=T)




#################### DRAW FINAL POINTS ################# 
##ALL###
nclus=8
colrfin <- idcludat$Cluster
colrfin<- get_colors(colrfin, group.col = brewer.pal(nclus, "Set2"))

####Figure 3.2A
lines3d(xyz[,1],xyz[,2],xyz[,3],col="black", xlim=rg,ylim=rg,zlim=rg,phi =30,theta =30, box=TRUE)
for (i in 1:dim(td)[1]) {
  points3d(tds[i,1],tds[i,2],tds[i,3],col=colrfin[i], add=T, pch=1, size=5)  # end point
}
axes3d(edges = "bbox", labels = TRUE, tick = TRUE, nticks = 5, box = F, expand = 1.5, xlab="x",ylab="y",zlab="z")
rgl.postscript("../", fmt = "pdf", drawText = TRUE )


#For plotting the final no colour, swap out above - Figure 3.1B
#points3d(tds[i,1],tds[i,2],tds[i,3],col="Black", add=T, pch=1, size=5)



#####INITIAL##############  - Figure 3.1A
decorate3d(xyz.coords(ctre[1,1],ctre[1,2],ctre[1,3]), expand = 1.5)
axes3d(edges = "bbox", labels = TRUE, tick = TRUE, nticks = 5, box = FALSE, expand = 1.25, xlab="x",ylab="y",zlab="z")
lines3d(xyz[,1],xyz[,2],xyz[,3],col="black", xlim=rg,ylim=rg,zlim=rg,phi =30,theta =30, box=T)

for (i in 1:dim(td)[1]) {
  points3d(td[i,1,1],td[i,1,2],td[i,1,3],col="black",add=T,pch=20, size=2)  # initial point
}
view3d(userMatrix = rglpos)
rgl.postscript("../", fmt = "pdf", drawText = TRUE )





#Subsetting 

lines3d(xyz[,1],xyz[,2],xyz[,3],col="black", xlim=rg,ylim=rg,zlim=rg,phi =30,theta =30, box=T)


##plot the initials by colour of final group - you can easily isolate a specific group
colll<- get_colors(c(1:8), group.col = brewer.pal(nclus, "Set2"))
lines3d(xyz[,1],xyz[,2],xyz[,3],col="black", xlim=rg,ylim=rg,zlim=rg,phi =30,theta =30, box=T)
points3d(cluin1[,2],cluin1[,3],cluin1[,4],col=colll[1], add=T, pch=19, cex=0.5, size=3)
points3d(cluin2[,2],cluin2[,3],cluin2[,4],col=colll[2], add=T, pch=19, cex=0.5, size=3)
points3d(cluin3[,2],cluin3[,3],cluin3[,4],col=colll[3], add=T, pch=19, cex=0.5, size=3)
points3d(cluin4[,2],cluin4[,3],cluin4[,4],col=colll[4], add=T, pch=19, cex=0.5, size=3)
points3d(cluin5[,2],cluin5[,3],cluin5[,4],col=colll[5], add=T, pch=19, cex=0.5, size=3)
points3d(cluin6[,2],cluin6[,3],cluin6[,4],col=colll[6], add=T, pch=19, cex=0.5, size=3)
points3d(cluin7[,2],cluin7[,3],cluin7[,4],col=colll[7], add=T, pch=19, cex=0.5, size=3)
points3d(cluin7[,2],cluin7[,3],cluin7[,4],col=colll[7], add=T, pch=19, cex=0.5, size=3)
axes3d(edges = "bbox", labels = TRUE, tick = TRUE, nticks = 5, box = FALSE, expand = 1.25, xlab="x",ylab="y",zlab="z")
rgl.postscript("../.pdf", fmt = "pdf", drawText = TRUE )



#GROUP BY FINAL CLUSTER

#INITAL SUBSET by FINAL
cluin1
cluin2
cluin3
cluin4
cluin5 ####NADP binding site
cluin6
cluin7
cluin8

#####SAMPLE TRAJ
collis
#Colour pal assigned to num
IDC1<- rownames(cluin1)
IDC2<-rownames(cluin2)
IDC3<-rownames(cluin3)
IDC4<- rownames(cluin4)
IDC5<-rownames(cluin5)####NADP binding site
IDC6<-rownames(cluin6)
IDC7<-rownames(cluin7)
IDC8<-rownames(cluin8)

####Novel trajectories were avoided for a simplier plot
#seleID1 1:15
#seleID2 1:15
#seleIDC3 1:4, 6:15
#seleIDC4 c(2,3,4,6,7,9,10,13,14,15,16,17,20,21,22)
#IDC5[c(1,3,4,5,6,7,9,10,11,13,14,15,16,17)
#IDC6[c(1,2,3,4,5,6,7,8,9,11,12,13,16,17)
#IDC7[c(1:15)]
#IDC8c(1,2,3,4,5,6,7,8,9,10,12,13,14,15,16)

##plotting of Appendix Figure 7.2
lines3d(xyz[,1],xyz[,2],xyz[,3],col="black", xlim=rg,ylim=rg,zlim=rg,phi =30,theta =30, box=T)
axes3d(edges = "bbox", labels = TRUE, tick = TRUE, nticks = 5, box = FALSE, expand=1.25, xlab="x",ylab="y",zlab="z")
for (i in as.numeric(IDC8[c(1,2,3,4,5,6,7,8,9,10,12,13,14,15,16)])) {
  points3d(td[i,1,1],td[i,1,2],td[i,1,3],col=collis[8], add=T, pch=1, cex=0.5)  # Initial point
  points3d(tds[i,1],tds[i,2],tds[i,3],col=collis[8], add=T, pch=1, size=5)
  lines3d(td[i,,1],td[i,,2],td[i,,3],add=T,col=collis[8])
}
rgl.postscript("../", fmt = "pdf", drawText = TRUE )
clear3d()












