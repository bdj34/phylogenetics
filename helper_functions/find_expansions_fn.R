find_expansions_fn <- function(tree, min_tips = 5, min_age = 5) {
  
  # Take time based ultrametric tree with edge lengths in years, 
  # find expansions with at least min_tips samples that happen after min_age years
  
  library(castor)
  library(ape)
  library(phangorn)
  
  #tree <- all_trees[[6]]$ultratree
  #min_age <- 5
  #min_tips <- 5
  
  #tree$tip.label
  
  root_node <- find_root(tree)
  
  dist_nodes_mat <- dist.nodes(tree)
  dist_from_root <- dist_nodes_mat[,root_node]
  
  # Find nodes happening after min_age, then check descendants of those
  recent_nodes <- which(dist_from_root > min_age)
  descendants_list <- Descendants(tree)
  all_descendants_list <- allDescendants(tree)
  names(descendants_list) <- as.character(seq(1, length(descendants_list)))
  names(all_descendants_list) <- as.character(seq(1, length(all_descendants_list)))
  descendants_of_recent_nodes <- descendants_list[as.character(recent_nodes)]
  all_descendants_of_recent_nodes <- all_descendants_list[as.character(recent_nodes)]
  
  # Find which have more than min_tips descendants
  n_descendants <- sapply(descendants_of_recent_nodes, length)
  nodes_keep_tips <- descendants_of_recent_nodes[n_descendants >= min_tips]
  nodes_keep_all_descendants <- all_descendants_of_recent_nodes[n_descendants >= min_tips]
  
  # Exclude those that are descendants of each other
  nodes_char <- names(nodes_keep_tips)
  descendants_of_selected_nodes <- c()
  for (i in 1:length(nodes_char)){
    descendants_of_selected_nodes <- c(descendants_of_selected_nodes, 
                                       unlist(nodes_keep_all_descendants[nodes_char[i]]))
  }
  final_nodes <- nodes_char[!nodes_char %in% as.character(descendants_of_selected_nodes)]
  
  # Make final list of trees
  tipsKeep <- nodes_keep_tips[final_nodes]
  subtree_list <- lapply(tipsKeep, keep.tip, phy = tree)
  names(subtree_list) <- names(tipsKeep)
  
  # Check that each node has distinct descendant tips
  stopifnot(length(unlist(tipsKeep)) == length(unique(unlist(tipsKeep))))
  
  return(subtree_list)

}
