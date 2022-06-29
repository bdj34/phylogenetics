init_growth_rate_vec <- function(birthRate_init_vec, nIters, 
                                                     popsize = 1e5, ndriver = 1, basefit = 0.2, 
                                                     nYears = 10, mutRate = 10000, nTips = 50) {
  
  out.df <- data.frame("lambda" = NULL, 
                       "method" = NULL,
                       "iteration" = NULL)
  
  for (birthRate_init in birthRate_init_vec) {
    
    lambda_NO_zeroes_vec <- c()
    lambda_zeroes_vec <- c()
    for (i in c(1:nIters)) {
      
      cfg_init = getDefaultConfig(target_pop_size  = 1e5, ndriver = 1, basefit = 0.2, rate = birthRate_init)
      sp_init = sim_pop(NULL,params=list(n_sim_days=365*10, b_stop_at_pop_size=1),cfg=cfg_init)
      
      popsize.df <- data.frame("time" = sp_init$timestamp, "popsize" = sp_init$pop.size)
      #ggplot(popsize.df, aes(x = time, y = log(popsize, base = 10))) + geom_point() + xlab("time (days)")
      #ggplot(popsize.df, aes(x = time, y = popsize)) + geom_point() + xlab("time (days)")
      
      ##Subsample tree and calculate growth rate
      sampledtree_init = get_subsampled_tree(sp_init, nTips) 
      sampledtree_init_mut = get_elapsed_time_tree(sampledtree_init, mutrateperdivision=0, backgroundrate = mutRate)
      #plot_tree(sampledtree_init, cex.label = 0.5)
      #t = plot_tree(sampledtree_init_mut, cex.label = 0.5)
      #node_labels(t,cex=0.5)
      
      # Run rtree fit? NO! Not working for now (not sure why)
      #agedf <- data.frame("tip.label" = sampledtree_init_mut$tip.label, "age" = c(0.00000000000000001, rep(max(sp_init$timestamp), 50)))
      #sampledtree_init_mut$agedf <- agedf
      #ultra_tree_init <- fit_tree(sampledtree_init_mut, switch_nodes = c(), model = "poisson_tree", niter = 1000)
      #plot_tree(ultra_tree_init$ultratree, cex.label = 0)
      
      # Subset to exclude zeroes outgroup and compare which is better
      keepTips <- sampledtree_init_mut$tip.label[sampledtree_init_mut$tip.label != "s1"]
      no_zeroes_mut_tree <- keep.tip(sampledtree_init_mut, keepTips)
      #plot_tree(no_zeroes_mut_tree)
      lambda_mut_init_no_zeroes <- ultra_tree_lambda_fn(no_zeroes_mut_tree, withZeroes = F)
      #print(paste0("Mutation tree lambda NO zeroes: ", lambda_mut_init_no_zeroes[1]*mutRate))
      
      lambda_mut_init <- ultra_tree_lambda_fn(sampledtree_init_mut, withZeroes = T)
      #print(paste0("Mutation tree lambda with zeroes: ", lambda_mut_init[1]*mutRate))
      
      #lambda_ultra_init <- ultra_tree_lambda_fn(ultra_tree_init$ultratree)
      #print(paste0("Ultrametric tree lambda: ", lambda_ultra_init[1]))
      #print(paste0("Actual lambda: ", birthRate_init))
      
      # Get growth rate (adjust for mutation rate) and add to vec
      lambda_NO_zeroes_vec <- c(lambda_NO_zeroes_vec, lambda_mut_init_no_zeroes[1]*mutRate)
      lambda_zeroes_vec <- c(lambda_zeroes_vec, lambda_mut_init[1]*mutRate)
      
    }
    
    add.df <- data.frame("lambda" = c(lambda_zeroes_vec, lambda_NO_zeroes_vec), 
                         "method" = c(rep("with_zeroes", nIters), rep("no_zeroes_group", nIters)),
                         "iteration" = rep(c(1:nIters), 2), "actual_lambda" = rep(birthRate_init, nIters*2))
    out.df <- rbind(out.df, add.df)
  }
  
  
  return(out.df)
}
