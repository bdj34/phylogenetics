site_freq_ultra_trees_fn <- function(tree) {
  
  # Note this requires a tree with a root that has two descendants (ie. doesn't go back to birth)
  
  descendant_list <- allDescendants(tree)
  descendant_df <- data.frame("Node" = (length(tree$tip.label)+2):max(tree$edge), "Parent" = NA,
                              "Edge_length" = NA, "n_cells" = NA)
  
  # Then find parent and edge length corresponding to each node
  dist_node_mat <- dist.nodes(tree)
  for (k in descendant_df$Node) {
    descendant_df$n_cells[descendant_df$Node == k] <- sum(descendant_list[[k]] < length(tree$tip.label)+1)
    descendant_df$Parent[descendant_df$Node == k] <- tree$edge[which(tree$edge[,2] == k ),1]
    descendant_df$Edge_length[descendant_df$Node == k] <- dist_node_mat[descendant_df$Parent[descendant_df$Node == k], k]
    # Check edge length to make sure it's the same using dist.nodes or tree itself
    stopifnot(descendant_df$Edge_length[descendant_df$Node == k] == tree$edge.length[which(tree$edge[,2] == k )])
  }
  
  # Finally, sum up and record the edge lengths which are shared by the same number of cells
  site_freq <- data.frame("n_cells" = unique(descendant_df$n_cells), "freq" = NA)
  for (k in site_freq$n_cells) {
    site_freq$freq[site_freq$n_cells == k] <- sum(descendant_df$Edge_length[descendant_df$n_cells == k])
  }
  
  return(site_freq)
  
}