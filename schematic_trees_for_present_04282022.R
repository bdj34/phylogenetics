library(ape)
library(phangorn)
library(phytools)
library(geiger)
library("Biostrings")
library("ggplot2")
library("ggtree")
library(treeio)
library(gprofiler2)
library(treemut)
library(adephylo)
library(castor)

rm(list=ls())

######################## Simulate a birth death process ########################

n_selected <- 20; tumorAge <- 40
my_tree <- rlineage(birth = .095, death = 0.02, Tmax = tumorAge)
my_tree_extant <- rbdtree(birth = .09, death = 0, Tmax = tumorAge)

color <- rep("black", length(my_tree$edge[,1])+1)
linetype <- rep("solid", length(my_tree$edge[,1])+1)
size <- rep(.5, length(my_tree$edge[,1])+1)

nodePathAllSelected <- c()
nodePathList <- nodepath(my_tree)

#Check tree 
plot_tree(my_tree_extant, cex.label = F, default_edge_color = "black")
ggtree(my_tree, color = color, linetype = linetype, size =size)+
  layout_dendrogram()

tipsKeep <- sample(c(1:length(my_tree_extant$tip.label)), n_selected)
subset_tree <- get_subtree_with_tips(my_tree_extant, only_tips = tipsKeep, force_keep_root = T)
subtree <- subset_tree$subtree

plot_tree(subtree, cex.label = 0, default_edge_color = "black", lwd = 2)
ggtree(subtree, color = "black", linetype = "solid", size = 0.5)+
  layout_dendrogram()

color <- c()
for (i in c(1:(2*n_selected-2))) {
  if (subtree$edge[i,2] < length(subtree$tip.label) + 1) {
    color <- c(color, "#EE442F") 
  } else {
    color <- c(color, "#63ACBE") 
  }
}
#color <- c(rep("#63ACBE", n_selected), rep("#EE442F", n_selected-1))
#ggtree(subtree, color = color, linetype = "solid", size = 0.5)+
#  layout_dendrogram()
plot_tree(subtree, cex.label = 0, default_edge_color = color, lwd = 2)









#Choose random node (make sure it's alive)
#for(i in 1:n_selected){
#  while(check <- F){
#    rand_node <- sample(1:length(nodePathList), 1, replace=F)
#    nodePath <- nodePathList[[rand_node]]
#    sum(my_tree$edge.length[nodePath])
#  }
#}
#nodePathList[[sample(1:length(nodePathList), 1, replace=F)]]

#Choose selected nodes and keep vec of their paths
for (i in 1:n_selected){
  nodePathSelected <- nodePathList[[sample(1:length(nodePathList), 1, replace=F)]]
  nodePathSelected <- nodePathSelected[nodePathSelected < 2*n_selected-1]
  
  nodePathAllSelected <- c(nodePathAllSelected, nodePathSelected)
}

sumCheck <- 0
for (i in 1:length(nodePathSelected)){
  if(i %% 2 != 0){
    next
  }else{
    sumCheck <- sumCheck + my_tree$edge.length[nodePathSelected]
  }
}

linetype[nodePathAllSelected] <- "solid"
size[nodePathAllSelected] <- 1
color[nodePathAllSelected] <- "black"

#Check tree 
ggtree(my_tree, color = color, linetype = linetype, size =size)+
  layout_dendrogram()


plot_tree(my_tree)


color <- rep("grey85", length(my_tree$edge)/2+1)
linetype <- rep("dotted", length(my_tree$edge)/2+1)
size <- rep(.5, length(my_tree$edge)/2+1)
nodePathAllSelected <- c()
nodePathList <- nodepath(my_tree)

#Choose selected nodes and keep vec of their paths
for (i in 1:n_selected){
  nodePathSelected <- nodePathList[[rand_node_extant[i]]]
  nodePathSelected <- nodePathSelected[nodePathSelected < 2*n-1]
  
  nodePathAllSelected <- c(nodePathAllSelected, nodePathSelected)
}

linetype[nodePathAllSelected] <- "solid"
size[nodePathAllSelected] <- 1

#Plot
ggtree(my_tree, branch.length="none", color = color, linetype = linetype, size =size)+
  layout_dendrogram()
ggtree(my_tree, branch.length="none", color = color, linetype = linetype, size =size,
       layout = "circular")

#Count Non-singletons and singletons and color accordingly
nodePathNonSingletons <- table(nodePathAllSelected)>1
nodeNamesNonSingletons <- as.numeric(names(nodePathNonSingletons))[as.vector(nodePathNonSingletons)]
nodeNamesSingletons <- as.numeric(names(nodePathNonSingletons))[!as.vector(nodePathNonSingletons)]
color[nodeNamesNonSingletons] <- "green"
color[nodeNamesSingletons] <- "red"

ggtree(my_tree, branch.length="none", color = color, linetype = linetype, size =size)+
  layout_dendrogram()

