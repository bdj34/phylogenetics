# Compare our lambda estimates to those of the phylofit using only the coal times

library("rsimpop")
library(ape)

rm(list = ls())

source("~/rotation_fall2021/phylofit_fns.R")#
source("~/rotation_fall2021/find_expansions_fn.R")
source("~/rotation_fall2021/ultra_tree_lambda_fn.R")

dir.create("~/rotation_fall2021/mitchell", showWarnings = F)

min_clade_size <- 10

# Read in tree
all_trees <- readRDS("~/rotation_fall2021/ultra_trees_mitchell_07012022.rds")#
#dat=read.table("~/rotation_fall2021/selection_coeff_all_mitchell.csv", header = T, sep = ",")#
#dat=dat %>% filter(number_samples >= min_clade_size)


for (ID in names(all_trees)){
  
  # Skip if younger than 5
  if (all_trees[[ID]]$agedf$age[1] < 5){
    next
  }
  
  ultra_tree <- all_trees[[ID]]$ultratree
  tree <- plot_tree(ultra_tree, b_do_not_plot = TRUE)
  #plot_tree(tree)
  numcolonies=length(tree$tip.label)
  
  expansions <- find_expansions_fn(ultra_tree, min_tips = min_clade_size)
  if (length(expansions) == 0) {
    next
  }
  expanded_nodes <- as.numeric(names(expansions))
  
  nh <- nodeHeights(tree)
  Age <- ultra_tree$agedf$age[1]
  
  out=lapply(1:length(expanded_nodes), function(i) {
    #i <- 1
    ht=nh[which(tree$edge[,2]==expanded_nodes[i]),2]
    t1=drop.tip(tree,setdiff(tree$tip.label,get_samples_in_clade(expanded_nodes[i],tree)))
    
    #t1_plot <- plot_tree(t1)
    #ggtree(tree)+geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + layout_dendrogram()
    
    our_lambda <- ultra_tree_lambda_fn(t1, withZeroes = F)[1]#
    
    t1=add_outgroup(t1,num.shared.var = ht)
    t1$coords=NULL
    plot_tree(t1,cex.label = 0,mar=par("mar"))
    
    ## Fit without taking into account final aberrant cell fraction
    res1=fit_clade(tree, expanded_nodes[i], nmutcolony = -1, nwtcolony =-1, maxt=1.0*Age,
                   minLN = 4, maxLN = 6, nchain = 6)
    #plot_res(res = res1,Age,b_add_extra = TRUE);title(sprintf('%s:No ACF Target',dat$variant_ID[i]))
    
    ## Fit taking into account final aberrant cell fraction
    #res2=fit_clade(tree,dat$node[i],nmutcolony = dat$number_samples[i],nwtcolony = numcolonies-dat$number_samples[1],maxt=1.5*Age)
    #plot_res(res = res2,Age,b_add_extra = TRUE);title(sprintf('%s:With ACF Target',dat$variant_ID[i]))
    list(res1=res1, our_lambda = our_lambda, subtree = t1)
  })
  
  saveRDS(out, paste0("~/rotation_fall2021/mitchell/out_lists_compare_phylofit_no_ACF_our_lambda_pt_", ID, "_07062022.rds"))
  
  df <- as.data.frame(t(as.data.frame(lapply(1:length(out), function(i){
    conf <- 0.95
    ptile=c(0.5*(1-conf),0.5,0.5*(1+conf))
    quantile(out[[i]]$res$posterior$S,ptile)
  }))))
  
  df$our_lambda <- c(unlist(lapply(1:length(out), function(i){return(out[[i]]$our_lambda)})))
  row.names(df) <- expanded_nodes
  df$node <- factor(expanded_nodes)
  colnames(df) <- c("lb", "median", "ub", "our_lambda", "node")
  df$N_tips <- c(unlist(lapply(1:length(out), function(i){return(length(out[[i]]$subtree$tip.label))})))
  
  pdf(paste0("~/rotation_fall2021/mitchell/plot_compare_phylofit_no_ACF_our_lambda_pt_", ID, "_07062022.pdf"))
  print(ggplot(df) + geom_point(aes(x = node, y = median), color = "blue") +
          geom_point(aes(x = node, y = our_lambda), color = "red") + 
          geom_errorbar(aes(x = node, ymin = lb, ymax = ub), color = "blue") +
          geom_text(aes(x = node, y = ub + 0.1, label = paste0("n = ", N_tips)))
  )
  dev.off()
  
  for ( i in 1:length(out)){
    pdf(paste0("~/rotation_fall2021/mitchell/subtree_pt_", ID, "_node_", expanded_nodes[i], "_07062022.pdf"))
    print(plot_tree(out[[i]]$subtree))
    dev.off()
  }
  
  
}
