mut_tree_lambda_fn <- function(subtree, age_post_conception) {
  
  # Only works for trees with the original zero root included
  
  #(nu*n)/(M_i) 
  
  # Estimate nu using the total number of mutations and the time since birth
  nu_vec <- distRoot(subtree)/age_post_conception
  
  # Get mutation rate estimates
  nu_mean <- mean(nu_vec)
  nu_sd <- sd(nu_vec)
  nu_max <- max(nu_vec)
  nu_min <- min(nu_vec)
  
  # Use subset tree to get the edge lengths corresponding to the number of descendants
  # First get list of descendents from each internal node
  descendant_list <- allDescendants(subtree)
  descendant_df <- data.frame("Node" = (length(subtree$tip.label)+3):max(subtree$edge), "Parent" = NA,
                              "Edge_length" = NA, "n_cells" = NA)
  
  # Then find parent and edge length corresponding to each node
  dist_node_mat <- dist.nodes(subtree)
  for (k in descendant_df$Node) {
    descendant_df$n_cells[descendant_df$Node == k] <- sum(descendant_list[[k]] < length(subtree$tip.label)+1)
    descendant_df$Parent[descendant_df$Node == k] <- subtree$edge[which(subtree$edge[,2] == k ),1]
    descendant_df$Edge_length[descendant_df$Node == k] <- dist_node_mat[descendant_df$Parent[descendant_df$Node == k], k]
    # Check edge length to make sure it's the same using dist.nodes or tree itself
    stopifnot(descendant_df$Edge_length[descendant_df$Node == k] == subtree$edge.length[which(subtree$edge[,2] == k )])
  }
  
  #write.csv(descendant_df, paste0(outDir, "/descendant_info_", 
  #                                patient_name, "_age_", round(age,2), "_driver_", gsub(":", "_", driver), ".csv"))
  
  # Finally, sum up and record the edge lengths which are shared by the same number of cells
  site_freq <- data.frame("n_cells" = unique(descendant_df$n_cells), "freq" = NA)
  for (k in site_freq$n_cells) {
    site_freq$freq[site_freq$n_cells == k] <- sum(descendant_df$Edge_length[descendant_df$n_cells == k])
  }
  
  # The sum of edge lengths in site_freq dataframe should equal the total internal branch lengths
  site_freq_internal <- sum(site_freq$freq)
  
  growthRate_mean <- nu_mean*length(subtree$tip.label)/site_freq_internal
  growthRate_mean_minus_sd <- (nu_mean-nu_sd)*length(subtree$tip.label)/site_freq_internal
  growthRate_mean_plus_sd <- (nu_mean+nu_sd)*length(subtree$tip.label)/site_freq_internal
  growthRate_min <- nu_min*length(subtree$tip.label)/site_freq_internal
  growthRate_max <- nu_max*length(subtree$tip.label)/site_freq_internal
  
  out_vec <- c(growthRate_mean, growthRate_mean_minus_sd, growthRate_mean_plus_sd, 
    growthRate_min, growthRate_max, site_freq_internal, nu_mean, nu_sd, nu_min, nu_max)
  
  names(out_vec) <- c("mean_growth_rate", "mean_growth_rate_minus_sd", "mean_growth_rate_plus_sd",
  "min_growth_rate", "max_growth_rate", "internal_lengths", "mean_mutRate",
  "sd_mutRate", "min_mutRate", "max_mutRate")
  
  return(out_vec)
}
