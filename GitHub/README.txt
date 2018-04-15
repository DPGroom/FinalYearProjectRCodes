README -

All code is included for analysis outlined within 2.2 as well as codes for the plotting the figures. 
The Codes remain extensively commented to guide the reader through the methods of analysis. 
Some code is co-dependent on products of one another; those of which have been highlighted below.

– “Simple 3D plotting.R” is self-contained and was used in the results outlined in 3.1-3.2. 
-	Initially the data is imported, and the initial and final frames are separated within matrixes 
	after removing non-binding trajectories. These are then plotted representing Figure 3.1. The Final 
	Coordinates are then grouped and then used to group the Initial coordinates which are then plotted 
	within Figure 3.2.
-	Following this trajectories are then also isolated and representative trajectories plotted for Appendix 
	Figure 7.2.

– “Network Analysis.R” is the master file that extensively covers the majority of this project and its products form the basis for this studies later analysis. 
-	Initially the data is imported and a KDE is done to identify regions of high spatial density. 
-	This is then followed by resolving individual Cluster-points and a kmeans cluster analysis conducted
	to produce basins of attraction.
-	Principal Coordinate Analysis is then conducted to spatially group the basins into diffusive binding sites. 
-	Trajectories are then assigned to basins and mobile trajectories (those that migrate to multiple basins) 
	are sample
-	The Network Edge List is constructed from these mobile trajectories 
-	A threshold of 0.005 is placed upon edges, those that fall below are disregarded to map common diffusion
	and remove edges for a more insightful analysis
-	The Steady State Analysis conducted in Figure 3.6
-	The network is plotted initially in 2D for Figure 3.4 then subsequently 3D on the protein backbone 4.1.C

– “Graph Theory.R” is dependent on the output files produced within Network Analysis.R and forms the basis for
	Methods 2.2.5 and Figures 3.4.C and 3.5.
-	The data is initially loaded in from outputs of Network Analysis.R
-	Subsets the directional In/Out-Degrees
-	Produces Figure 3.5
-	The transition matrix is used to plot the heatmap present within Figure 3.4.C for a more pleasing 
	visualisation


– “RMSD.R” is also dependent on output files for Network Analysis.R for plotting Appendix Figure 7.5
-	RMSD data and previous cluster data from Network Analysis.R is imported
-	The Average RMSD for trajectories that visit a cluster is calculated per frame and then used to plot 
	Appendix Figure 7.5
-	A sample trajectory Appendix Figure 7.4 values are also taken and plotted 
