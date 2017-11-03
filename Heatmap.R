#Setup
rm(list=ls())

#Install libraries
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = TRUE) }
if (!require("Cairo")) { install.packages("Cairo", dependencies = TRUE) }
if (!require("RColorBrewer")) { install.packages("RColorBrewer", dependencies = TRUE) }
if (!require("ape")) { install.packages("ape", dependencies = TRUE) }
if (!require("gplots")) { install.packages("gplots", dependencies = TRUE) }

# creates own color palette from red to green
my_palette <- colorRampPalette(c( "red", "orangered", "orange", "yellow"))(n = 49)



####################### Loop ###############################

setwd("~/LSE_pipeline/Stage3/Spada/Rfigs")
path = "~/LSE_pipeline/Stage3/Spada/Rfigs"
file.names <- dir(path, pattern =".dnd")
file.names2<- dir(path,  pattern =".tbl2")

for(i in 1:23){
  #i=2
  #Read tree file
  tree <- read.tree(file.names[i])

  #Transform branch lengths of zero
  rep_tree_r <- root(tree, resolve.root = T, node=length(tree$tip.label)+4)
  rep_tree_r$edge.length[which(rep_tree_r$edge.length <= 0)] <- 0.00001
  rep_tree_um <- chronopl(rep_tree_r, lambda = 0.1, tol = 0)
  rep_tree_d <- as.dendrogram(as.hclust.phylo(rep_tree_um))
  #Get info on label order
  clade_order <- order.dendrogram(rep_tree_d)
  clade_name <- labels(rep_tree_d)
  clade_position <- data.frame(clade_name, clade_order)
  clade_position <- clade_position[order(clade_position$clade_order),]
  clade_order
  #Read in data and transform it into matrix format
  GroupMat <- read.delim(file.names2[i])
  rnames <- GroupMat[,1]                            # assign labels in column 1 to "rnames"
  mat_data <- data.matrix(GroupMat[,2:ncol(GroupMat)])  # transform column 2-12 into a matrix
  rownames(mat_data) <- rnames                  # assign row names 
  #Force row order to match the order of leafs in rep_tree_d
  new_order <- match(clade_position$clade_name, row.names(mat_data))
  mat_data <- mat_data[new_order,]
  #Reorder table to agree with SeqID order on the tree
  mat_data[which(mat_data < 0.01)] <- 0.01
  mat_data <- log(mat_data)
  #Customize and plot the heat map
  specificity <- as.vector(mat_data[,ncol(mat_data)])
  cc <- ifelse(specificity == 0, "darkorchid4", "white")
  my_plot1 <- function(data){heatmap(data[,1:ncol(data)-1], na.rm=TRUE, Colv = NA,  
                                     Rowv = rep_tree_d, margins = c(0,13), 
                                     cexRow=0.5, cexCol=1.2,
                                     scale = "none", col=my_palette, RowSideColors=cc
  )}

  my_plot1(mat_data)
 


#Save image
my_plot1 <- function(data){heatmap(data[,1:ncol(data)-1], na.rm=TRUE, Colv = NA,  
                                   Rowv = rep_tree_d, margins = c(13,130), 
                                   cexRow = 2 + 2.8/log10(dim(data)[1L]), cexCol = 2 + 1.2/log10(dim(data)[2L]),
                                   scale = "none", col=my_palette, RowSideColors=cc )}
png(filename=paste(file.names[i], ".png", sep=""), 
      type="cairo",
      units="in", 
      width=150, 
      height=50, 
      pointsize=15, 
      res=96)
  my_plot1(mat_data)
  dev.off()

}






