get_coal_times_from_ultra_fn <- function(tree) {
  
  # Note: This only works for the trees which remove the original zeros root
  
  # Get the internal nodes of the tree
  internal_nodes <- c((length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode))
  
  # Get the distance matrix
  dist_node_mat <- dist.nodes(tree)
  
  ## then get the edge lengths for internal nodes (coal_times)
  coal_times <- dist_node_mat[1,(length(tree$tip.label)+1)] - dist_node_mat[internal_nodes,(length(tree$tip.label)+1)]
  
  if (length(coal_times) != n-1) {
    # Find multifurcating node
    multi_node <- as.numeric(names(table(tree$edge[,1]))[which(table(tree$edge[,1]) == 3)])
    print(multi_node)
    internal_nodes_v2 <- c(multi_node, internal_nodes)
    coal_times <- dist_node_mat[1,(length(tree$tip.label)+1)] - dist_node_mat[internal_nodes_v2,(length(tree$tip.label)+1)]
    stopifnot(length(coal_times) == n-1)
  }
  
  return(coal_times)
}

coal_times_list_to_mat <- function(coal_times_list) {
  
  n.obs <- sapply(coal_times_list, length)
  seq.max <- seq_len(max(n.obs))
  mat <- t(sapply(coal_times_list, "[", i = seq.max))
  colnames(mat) <- paste0("Coal_time", c(1:dim(mat)[2]))
  return(mat)
}



