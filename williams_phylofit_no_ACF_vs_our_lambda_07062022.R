# Compare our lambda estimates to those of the phylofit using only the coal times
# Get growth rates from VCF trees before and after ultrametric conversion
# No differing mutation rates at birth or with drivers
library(ape)
library(phangorn)
library(treemut)
library(rtreefit)
library(ggtree)
library(ggplot2)
library(vcfR)
library(readr)
library(castor)
library(adephylo)

rm(list = ls())

source("~/rotation_fall2021/phylofit_fns.R")#
source("~/rotation_fall2021/find_expansions_fn.R")
source("~/rotation_fall2021/ultra_tree_lambda_fn.R")
source("~/rotation_fall2021/get_coal_times_from_ultra_fn.R")

dir.create("~/rotation_fall2021/williams", showWarnings = F)

date <- "07062022"

# Set threshold number of tips at which we will calculate clade age (default = 10)
threshold_n_tips <- 10

# Import information about samples including drivers, ultrametric trees from paper, etc.
all_ultra <- readRDS('~/Downloads/PDD_TELO.rds')

# Import clade info (from table in Ext. Data Fig 6B)
clade_info <- read.table("~/rotation_fall2021/clade_info.txt", header = T)
clade_info$lambda <- log(clade_info$lambda/100 + 1)
clade_info$lower_lambda <- log(clade_info$lower_lambda/100 + 1)
clade_info$upper_lambda <- log(clade_info$upper_lambda/100 + 1)

# Make df for saving
df <- data.frame("pid" = NA, "driver" = NA, "phylofit_growthRate" = NA, "phylofit_lb_growthRate" = NA, 
                 "phylofit_ub_growthRate" = NA, "our_growthRate_theirTree" = NA, 
                 "N_tips" = NA, "age_at_sample" = NA)

pids <- names(all_ultra)
count <- 1
coal_times_list_ultra_theirs <- list()

for (i in 1:length(pids)) {
  #i <- 1
  
  patientName <- pids[i]
  print(paste0("Running patient ", i, ": ", patientName))
  patientInfo <- all_ultra[[i]]
  
  # Get agedf and number of samples from unique (nonzero) timepoints
  agedf <- patientInfo$agedf
  
  #plot_tree(our_ultra_tree$ultratree)
  #ggtree(our_ultra_tree$ultratree) + layout_dendrogram() + geom_tiplab(angle = 90)
  
  #Fix notation in agedf metadata
  agedf$driver3 <- gsub(",", ":", agedf$driver3)
  
  # 5182 is poorly annotated, skip for now
  if (patientName == "PD5182") {
    next
  }
  
  # Go through each driver at each timepoint
  for (driver in unique(agedf$driver3)) {
    for (age in unique(agedf$age_at_sample_pcy)) {
      #driver <- unique(agedf$driver3)[2]
      #age <- unique(agedf$age_at_sample_pcy)[1]
      if (patientName == "PD9478" & driver == "JAK2"){
        driver <- "JAK2:DNMT3A"
        agedf$driver3[agedf$driver3=="JAK2"] <- "JAK2:DNMT3A"
      }
      
      n <- sum(agedf$driver3 == driver & agedf$age_at_sample_pcy == age)
      
      if (n >= threshold_n_tips & driver != "WT") {
        
        # Get indexes of the samples that have this age and driver combo
        rowsKeep <- which(agedf$driver3 == driver & agedf$age_at_sample_pcy == age)
        
        tipsKeep <- agedf$tip.label[rowsKeep]
        
        ########## Include coincidental CNV in PD5179 ??? ##########
        #if (patientName == "PD5179") {
        #  tipsKeep <- c(tipsKeep, "bo", "bs")
        #  n <- length(tipsKeep)
        #}
        
        # Get growth rate from their ultrametric tree (from theirTree)
        theirTree <- patientInfo$ultratree
        subtree_theirTree <- keep.tip(theirTree, tipsKeep)
        #ggtree(subtree_theirTree)+layout_dendrogram()+geom_nodepoint()
        #ggtree(theirTree) + geom_text2(aes(label=node), hjust=-.3) + layout_dendrogram()
        out_theirTree <- ultra_tree_lambda_fn(subtree_theirTree)
        growthRate_theirTree <- out_theirTree[1]
        internal_lengths_theirTree <- out_theirTree[2]
        # Get coalescence times from their ultrametric tree
        coal_times_ultra_theirs <- get_coal_times_from_ultra_fn(subtree_theirTree)
        
        driver_node <- getMRCA(theirTree, tipsKeep)
        res1=fit_clade(theirTree, driver_node, nmutcolony = -1, nwtcolony =-1, maxt=1.0*age,
                       minLN = 4, maxLN = 6, nchain = 6)
        
        conf <- 0.95
        ptile=c(0.5*(1-conf),0.5,0.5*(1+conf))
        their_result <- quantile(res1$posterior$S,ptile)
        
        # Determine if driver is in the inputted clade_info table (need their growth rate)
        bool <- F
        for (l in 1:length(clade_info$driver)) {
          bool <- all(unlist(strsplit(driver, ":")) %in% unlist(strsplit(clade_info$driver[l], ":"))) &
            all(unlist(strsplit(clade_info$driver[l], ":")) %in% unlist(strsplit(driver, ":"))) &
            clade_info$pid[l] == patientName
          if (bool) {
            row <- l
            break
          }            
        }
        
        if (bool) {
          
          df[count,] <- c(clade_info[row,c(1:2)], their_result[2], their_result[1], 
                          their_result[3], growthRate_theirTree, n, age)
          
          coal_times_list_ultra_theirs <- append(coal_times_list_ultra_theirs, list(coal_times_ultra_theirs))
          
          count <- count + 1
          #print(paste0("Our tumor age: ", ourAge, "   their tumor age: ", theirAge, "   their tumor lambda: ", lambda, "   fit lambda: ", our_lambda_their_age))
          
        }else {
          pid_clades <- clade_info[clade_info$pid == patientName,]
          print(paste0("In patient ", patientName, ", ", driver, " not in ", pid_clades$driver))
        }
      }
    }
  }
}

write.csv(df, paste0("~/rotation_fall2021/williams_clade_estimates_threshold_n_", threshold_n_tips,"_phylofit_no_ACF", date, ".csv"))

plot.df <- df
plot.df$driver_and_pid_and_age <- paste0(plot.df$pid, " ", plot.df$driver, " ", round(plot.df$age_at_sample, 2))
plot.df$driver_and_pid <- paste0(plot.df$pid, " ", plot.df$driver)

# Write coal_times to csv
names(coal_times_list_ultra_theirs) <- plot.df$driver_and_pid_and_age
coal_times_mat_ultra_theirs <- coal_times_list_to_mat(coal_times_list_ultra_theirs)
row.names(coal_times_mat_ultra_theirs) <- names(coal_times_list_ultra_theirs)

write.csv(coal_times_mat_ultra_theirs, paste0("~/rotation_fall2021/coal_times_from_their_trees_threshold_n_", threshold_n_tips, "_", date, ".csv"))

# Set plot ordering
plot.df$driver_and_pid <- factor(plot.df$driver_and_pid, 
                                 levels = unique(plot.df$driver_and_pid[order(plot.df$phylofit_growthRate)]))
plot.df$driver_and_pid_and_age <- factor(plot.df$driver_and_pid_and_age, 
                                         levels = unique(plot.df$driver_and_pid_and_age[order(plot.df$phylofit_growthRate)]))

ggplot(plot.df, aes(x = driver_and_pid_and_age)) + geom_point(aes(y = phylofit_growthRate), color = "blue") + 
  geom_errorbar(aes(ymin=phylofit_lb_growthRate, ymax=phylofit_ub_growthRate), width=.2, color = "blue") + 
  geom_point(aes(y = our_growthRate_theirTree), color = "red", size = 2, alpha = .8) + #, fill = "black", shape = 22) +
  geom_text(aes(label = paste0("n=", N_tips), y = 2.0))+
  ylab("Growth Rates (per year)") + #ylim(0, 2) + theme_bw() +
  ggtitle("Phylofit growth rate (blue) vs. our growth rates (red)") +
  theme(axis.text.x = element_text(face = "bold", angle = 50, vjust = 1, hjust=1, size = 10),
        panel.grid.minor = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 10))


ggplot(plot.df, aes(x = driver_and_pid_and_age)) + geom_point(aes(y = mutRate_rtreefit_whole_tree_mean), color = "blue") + 
  geom_errorbar(aes(ymin=mutRate_rtreefit_whole_tree_min, ymax=mutRate_rtreefit_whole_tree_max), width=.2, color = "blue") + 
  geom_point(aes(y = mutRate_VCF_clade_mean), color = "red", fill = "black", size = 2, alpha = .8) +
  geom_errorbar(aes(ymin=mutRate_VCF_clade_mean - mutRate_VCF_clade_sd, 
                    ymax=mutRate_VCF_clade_mean + mutRate_VCF_clade_sd), width=.2, color = "red")+
  geom_text(aes(label = paste0("n=", N_tips), y = 2.0))+
  ylab("Mutation Rates (per year)") + theme_bw() +
  ggtitle("Mutation Rates (mean +/- sd): Rtreefit rate for whole tree (blue) vs. mutation counting in clade (red square)") +
  theme(axis.text.x = element_text(face = "bold", angle = 70, vjust = 1, hjust=1), 
        panel.grid.minor = element_blank(), axis.title.x = element_blank())























library("rsimpop")
library(ape)

rm(list = ls())

source("~/rotation_fall2021/phylofit_fns.R")#
source("~/rotation_fall2021/find_expansions_fn.R")
source("~/rotation_fall2021/ultra_tree_lambda_fn.R")

dir.create("~/rotation_fall2021/williams", showWarnings = F)

min_clade_size <- 10

# Read in tree
all_trees <- readRDS('~/Downloads/PDD_TELO.rds')
#all_trees <- readRDS("~/rotation_fall2021/ultra_trees_williams_07012022.rds")#
#dat=read.table("~/rotation_fall2021/selection_coeff_all_williams.csv", header = T, sep = ",")#
#dat=dat %>% filter(number_samples >= min_clade_size)


for (ID in names(all_trees)){
  
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
  Age <- ultra_tree$agedf$age[2]
  
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
  
  saveRDS(out, paste0("~/rotation_fall2021/williams/out_lists_compare_phylofit_no_ACF_our_lambda_pt_", ID, "_07062022.rds"))
  
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
  
  pdf(paste0("~/rotation_fall2021/williams/plot_compare_phylofit_no_ACF_our_lambda_pt_", ID, "_07062022.pdf"))
  print(ggplot(df) + geom_point(aes(x = node, y = median), color = "blue") +
    geom_point(aes(x = node, y = our_lambda), color = "red") + 
    geom_errorbar(aes(x = node, ymin = lb, ymax = ub), color = "blue") +
    geom_text(aes(x = node, y = ub + 0.1, label = paste0("n = ", N_tips)))
  )
  dev.off()
  
  for ( i in 1:length(out)){
    pdf(paste0("~/rotation_fall2021/williams/subtree_pt_", ID, "_node_", expanded_nodes[i], "_07062022.pdf"))
    print(plot_tree(out[[i]]$subtree))
    dev.off()
  }
  
  
}

#ultra_tree <- all_trees[['KX008']]$ultratree
#expansions <- find_expansions_fn(ultra_tree, min_tips = min_clade_size)
#expanded_nodes <- as.numeric(names(expansions))
#ggtree(ultra_tree) + geom_cladelabel(node=expanded_nodes, label="Some random clade", color="red")+layout_dendrogram()



