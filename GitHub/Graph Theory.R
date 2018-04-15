library('RColorBrewer')
library(ggplot2)
library(ggfortify)
require(reshape2)
library(markovchain)
library(scales)


###Methods 2.2.5


##IMPORT THE NETRWORK
W_edgelist <- read.csv("Weighted Edge List - Rigid Vac.csv",sep=",")
W_edgelist$X<-NULL

#Processing of Markov without Network Analysis incase seperate runs
net2 <- graph_from_data_frame(W_edgelist)
ma_transitionProbs <- prop.table(as.matrix(get.adjacency(net2, attr="Probability")),1)
W_adjacency<-as.matrix(get.adjacency(net2, attr="Probability")) 
wammend2 <- as.data.frame(W_edgelist[W_edgelist[,3]>=0.005,,drop=F])
wammend2$Probability <- specify_decimal(wammend2$Probability,3)
###for plotting colours later on
wammendmark <- wammend2[wammend2$N1==wammend2$N2,,drop=F]
wammendsubset <- wammend2[wammend2$N1!=wammend2$N2,,drop=F]
net<- graph_from_data_frame(wammend2)


##QUICK LITTLE SUBSET OF COMMON DIFFUSION - seperates in and out degrees
nclus<-8
out<-NULL
inn<-NULL
for (i in 1:nclus){
  a <-wammend2[wammend2$N1==i,,drop=F]
  out <- rbind(out,a)
  b <-wammend2[wammend2$N2==i,,drop=F]
  inn <- rbind(inn,b)
}

###PROCESS IT with nice names
invout<-cbind(as.data.frame(table(t(out$N1))), as.data.frame(table(t(inn$N2)))$Freq, graph.strength(net))
colnames(invout)<- c("Cluster","Out","In","Total")
invout<-cbind(invout$Cluster,invout$Total,invout$In, invout$Out)
colnames(invout)<- c("Cluster","Total","In","Out")
invout<- as.data.frame(invout)
m.invout<- melt(invout, id.vars='Cluster')
colnames(m.invout)<- c("Cluster", "Edge", "value")

##CONNECTIVITY PLOT - Figure 3.5A
bar <- ggplot(m.invout, aes(Cluster, value))+   
  geom_bar(aes(fill = Edge), position = "dodge", stat="identity", color="black")+
  theme_classic()+
  scale_fill_brewer(palette = "Set2", direction = 1)+
  theme(legend.box.background = element_rect(),legend.box.margin = margin(4, 4, 4, 4))+
  ylab("Degree") + scale_y_continuous(breaks = c(2,4,6,8,10), expand = c(0, 0))+
  scale_x_continuous(breaks = c(1:nclus))+
  guides(color=guide_legend("Edge"))

##simple his - not used in report 
qplot(invout$In,geom="histogram",binwidth = 0.5,  xlab = "In-Degree")
qplot(invout$Out,geom="histogram",binwidth = 0.5,  xlab = "Out-Degree")

#extract legend
legend <- g_legend(bar) 
grid.draw(legend) 


###TOTAL  - Figure 3.5B
p1<-ggplot(SubsetTotal, aes(x=value))+
  geom_histogram(aes(y=..count../sum(..count..)), fill="#66c2b9", color="black", binwidth = 1)+
  theme_classic()+
  scale_fill_brewer(palette = "Set2", direction = 1)+
  ylab("") + scale_y_continuous(expand = c(0, 0))+
  xlab("Total Degree") + scale_x_continuous(breaks = c(0:10), expand = c(0, 0))


###IN - Figure 3.5B
p2<-ggplot(SubsetIn, aes(x=value))+
  geom_histogram(aes(y=..count../sum(..count..)), fill="#fc8d62", color="black", binwidth = 1)+
  theme_classic()+
  scale_fill_brewer(palette = "Set2", direction = 1)+
  ylab("Fraction of Nodes") + scale_y_continuous(expand = c(0, 0))+
  xlab("In-Degree") + scale_x_continuous(breaks = c(0:10), expand = c(0, 0))


##out - Figure 3.5B
p3<-ggplot(SubsetOut, aes(x=value))+
  geom_histogram(aes(y=..count../sum(..count..)), fill="#8da0cb", color="black", binwidth = 1)+
  theme_classic()+
  scale_fill_brewer(palette = "Set2", direction = 1)+
  ylab("") + scale_y_continuous(expand = c(0, 0))+
  xlab("Out-Degree") + scale_x_continuous(breaks = c(0:10), expand = c(0, 0))+ expand_limits(x = 6)

#DEGREE OF DISTRIBUTION his -  - Figure 3.5B
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size = "last"))


###HEAT MAP OF THE NETWORK -  - Figure 3.4C
hmp <- melt(ma_transitionProbs)
colnames(hmp) <- c("Var1","Var2","Probability")
jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "Blues"))
paletteSize <- 64
jBuPuPalette <- jBuPuFun(paletteSize)
ggplot(data=hmp, aes(x=Var1, Var2))+
  geom_tile(aes(fill=Probability),colour="White",size=1)+
  scale_fill_gradient2(low = jBuPuPalette[1],
                       mid = jBuPuPalette[paletteSize/2],
                       high = jBuPuPalette[paletteSize],
                       midpoint = (max(hmp$Probability) + min(hmp$Probability)) / 2,
                       name = "Probability")+
  theme_grey(base_size = 10) +
  theme(panel.background = element_blank(), 
        plot.title = element_text(size = 50, colour = "gray50"))+
  scale_x_continuous(breaks = c(1:nclus))+
  scale_y_continuous(breaks = c(1:nclus))+ xlab("Cluster")+ ylab("Cluster")

