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

date <- "06282022"

tree_dir <- paste0("~/rotation_fall2021/mutation_trees_mpBoot_treeMut/")
dir.create(tree_dir)
growthRate_dir <- paste0("~/rotation_fall2021/growthRates_", date, "/")
dir.create(growthRate_dir)
source("~/rotation_fall2021/build_mpBoot_tree_assign_edge_len_fn.R")
source("~/rotation_fall2021/ultra_tree_lambda_fn.R")
source("~/rotation_fall2021/mut_tree_lambda_fn.R")
source("~/rotation_fall2021/get_coal_times_from_ultra_fn.R")

# Set threshold number of tips at which we will calculate clade age (default = 10)
threshold_n_tips <- 5

# Import information about samples including drivers, ultrametric trees from paper, etc.
all_ultra <- readRDS('~/Downloads/PDD_TELO.rds')

# Import clade info (from table in Ext. Data Fig 6B)
clade_info <- read.table("~/rotation_fall2021/clade_info.txt", header = T)
clade_info$lambda <- log(clade_info$lambda/100 + 1)
clade_info$lower_lambda <- log(clade_info$lower_lambda/100 + 1)
clade_info$upper_lambda <- log(clade_info$upper_lambda/100 + 1)

# Make df for saving
df <- data.frame("pid" = NA, "driver" = NA, "their_growthRate" = NA, "their_lb_growthRate" = NA, 
                 "their_ub_growthRate" = NA, "our_growthRate_theirTree" = NA, 
                 "our_growthRate_VCF_ultra_Tree" = NA,  "our_growthRate_VCF_mut_Tree" = NA,
                 "our_growthRate_VCF_mut_Tree_minusSD" = NA, "our_growthRate_VCF_mut_Tree_plusSD" = NA,
                 "growthRate_mut_VCF_min" = NA, "growthRate_mut_VCF_max" = NA, 
                 "mutRate_VCF_clade_mean" = NA, "mutRate_VCF_clade_sd" = NA, "mutRate_VCF_clade_min" = NA,
                 "mutRate_VCF_clade_max" = NA, "mutRate_rtreefit_whole_tree_mean" = NA, 
                 "mutRate_rtreefit_whole_tree_median" = NA, "mutRate_rtreefit_whole_tree_sd" = NA,
                 "mutRate_rtreefit_whole_tree_min" = NA, "mutRate_rtreefit_whole_tree_max" = NA,
                 "N_tips" = NA, "internal_times_VCF_ultra" = NA, "internal_muts_VCF" = NA, 
                 "internal_times_theirTree" = NA, "age_at_sample" = NA)

pids <- names(all_ultra)
count <- 1
coal_times_list_ultra_VCF <- list()
coal_times_list_ultra_theirs <- list()

for (i in 1:length(pids)) {
  #i <- 8
  
  patientName <- pids[i]
  print(paste0("Running patient ", i, ": ", patientName))
  patientInfo <- all_ultra[[i]]
  
  vcf_file <- paste0("~/rotation_fall2021/vcfs_with_matched_sample_names/", patientName, ".rds")
  
  # Read VCF and organize
  vcf <- readRDS(vcf_file)
  data_mat_all <- vcf@gt[,2:dim(vcf@gt)[2]]
  mut_info <- vcf@fix[,"INFO"]
  row.names(data_mat_all) <- mut_info
  
  ref_nuc <- vcf@fix[,"REF"]
  alt_nuc <- vcf@fix[,"ALT"]

  tree_with_mutations <- build_mpBoot_tree_assign_edge_len_fn(vcf, tree_dir, patientName)
  
  #plot_tree(tree_with_mutations, cex.label = 0)
  #ggtree(tree_with_mutations) + layout_dendrogram() + geom_tiplab()

  # Get agedf and number of samples from unique (nonzero) timepoints
  agedf <- patientInfo$agedf
  n_samples <- length(unique(agedf$age_at_sample_pcy[agedf$age_at_sample_pcy > 1]))
  
  our_ultra_tree_file <- paste0(tree_dir, patientName, "ultrametric_tree_from_vcf.rds")
  
  if (file.exists(our_ultra_tree_file)) {
    our_ultra_tree <- readRDS(our_ultra_tree_file)
  } else {
    # Make tree ultrametric, assuming fixed mutation rate?
    age_df_for_tree_fit <- data.frame("tip.label" = agedf$tip.label, "age" = agedf$age_at_sample_pcy)
    tree_with_mutations$agedf <- age_df_for_tree_fit
    our_ultra_tree <- fit_tree(tree_with_mutations, switch_nodes = c(), model = "poisson_tree")
    saveRDS(our_ultra_tree, our_ultra_tree_file)
  }
  
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
        
        # Get growth rate from our mutation tree (from VCF)
        subset_tree_mut_VCF <- get_subtree_with_tips(tree_with_mutations, only_tips = tipsKeep, force_keep_root = T)
        subtree_mut_VCF <- subset_tree_mut_VCF$subtree
        #ggtree(subtree_mut_VCF) + layout_dendrogram() + geom_nodepoint()
        #ggtree(subsetTREE_MUTVCF) + layout_dendrogram() + geom_nodepoint()
        out_mut_VCF <- mut_tree_lambda_fn(subtree_mut_VCF, age)
        # Assign output to variables
        growthRate_mut_VCF <- out_mut_VCF[["mean_growth_rate"]]
        growthRate_mut_VCF_minus_sd <- out_mut_VCF[["mean_growth_rate_minus_sd"]]
        growthRate_mut_VCF_plus_sd <- out_mut_VCF[["mean_growth_rate_plus_sd"]]
        growthRate_mut_VCF_min <- out_mut_VCF[["min_growth_rate"]]
        growthRate_mut_VCF_max <- out_mut_VCF[["max_growth_rate"]]
        internal_muts_VCF <- out_mut_VCF[["internal_lengths"]]
        mutRate_VCF_mean <- out_mut_VCF[["mean_mutRate"]]
        mutRate_VCF_sd <- out_mut_VCF[["sd_mutRate"]]
        mutRate_VCF_min <- out_mut_VCF[["min_mutRate"]]
        mutRate_VCF_max <- out_mut_VCF[["max_mutRate"]]
        
        # Get growth rate from our ultrametric tree (from VCF)
        subtree_ultra_VCF <- keep.tip(our_ultra_tree$ultratree, tipsKeep)
        #ggtree(subtree_ultra_VCF)+geom_text2(aes(subset=!isTip, label=node), hjust=-.3)
        out_ultra_VCF <- ultra_tree_lambda_fn(subtree_ultra_VCF)
        growthRate_ultra_VCF <- out_ultra_VCF[1]
        internal_times_VCF <- out_ultra_VCF[2]
        # Get coalescence times from our ultrametric tree (from VCF)
        coal_times_ultra_VCF <- get_coal_times_from_ultra_fn(subtree_ultra_VCF)
        
        # Get growth rate from their ultrametric tree (from theirTree)
        subtree_theirTree <- keep.tip(patientInfo$ultratree, tipsKeep)
        #ggtree(subtree_theirTree)+layout_dendrogram()+geom_nodepoint()
        #ggtree(subtree_theirTree) + geom_text2(aes(label=node), hjust=-.3)
        out_theirTree <- ultra_tree_lambda_fn(subtree_theirTree)
        growthRate_theirTree <- out_theirTree[1]
        internal_lengths_theirTree <- out_theirTree[2]
        # Get coalescence times from their ultrametric tree
        coal_times_ultra_theirs <- get_coal_times_from_ultra_fn(subtree_theirTree)
        
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
          
          df[count,] <- c(clade_info[row,c(1:5)], growthRate_theirTree, growthRate_ultra_VCF, 
                          growthRate_mut_VCF, growthRate_mut_VCF_minus_sd, growthRate_mut_VCF_plus_sd, 
                          growthRate_mut_VCF_min, growthRate_mut_VCF_max, mutRate_VCF_mean,
                          mutRate_VCF_sd, mutRate_VCF_min, mutRate_VCF_max,
                          our_ultra_tree$lambda$mean, our_ultra_tree$lambda$median, 
                          our_ultra_tree$lambda$sd, our_ultra_tree$lambda$lb, our_ultra_tree$lambda$ub,
                          n, internal_times_VCF, internal_muts_VCF, internal_lengths_theirTree, age)
          
          coal_times_list_ultra_VCF <- append(coal_times_list_ultra_VCF, list(coal_times_ultra_VCF))
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

write.csv(df, paste0("~/rotation_fall2021/clade_estimates_threshold_n_5_", date, ".csv"))

plot.df <- df
plot.df$driver_and_pid_and_age <- paste0(plot.df$pid, " ", plot.df$driver, " ", round(plot.df$age_at_sample, 2))
plot.df$driver_and_pid <- paste0(plot.df$pid, " ", plot.df$driver)

# Write coal_times to csv
names(coal_times_list_ultra_VCF) <- plot.df$driver_and_pid_and_age
names(coal_times_list_ultra_theirs) <- plot.df$driver_and_pid_and_age
coal_times_mat_ultra_VCF <- coal_times_list_to_mat(coal_times_list_ultra_VCF)
coal_times_mat_ultra_theirs <- coal_times_list_to_mat(coal_times_list_ultra_theirs)
row.names(coal_times_mat_ultra_VCF) <- names(coal_times_list_ultra_VCF)
row.names(coal_times_mat_ultra_theirs) <- names(coal_times_list_ultra_theirs)

#write.csv(coal_times_mat_ultra_VCF, paste0("~/rotation_fall2021/coal_times_from_VCF_", date, ".csv"))
#write.csv(coal_times_mat_ultra_theirs, paste0("~/rotation_fall2021/coal_times_from_their_trees_", date, ".csv"))

# Set plot ordering
plot.df$driver_and_pid <- factor(plot.df$driver_and_pid, 
                                levels = unique(plot.df$driver_and_pid[order(plot.df$their_growthRate)]))
plot.df$driver_and_pid_and_age <- factor(plot.df$driver_and_pid_and_age, 
                                 levels = unique(plot.df$driver_and_pid_and_age[order(plot.df$their_growthRate)]))

ggplot(plot.df, aes(x = driver_and_pid_and_age)) + geom_point(aes(y = their_growthRate), color = "blue") + 
  geom_errorbar(aes(ymin=their_lb_growthRate, ymax=their_ub_growthRate), width=.2, color = "blue") + 
  geom_point(aes(y = our_growthRate_theirTree), color = "red", size = 2, alpha = .8) + #, fill = "black", shape = 22) +
  #geom_point(aes(y = our_growthRate_VCF_ultra_Tree), color = "red", fill = "black", shape = 24, size = 2, alpha = .8) + 
  #geom_point(aes(y = our_growthRate_VCF_mut_Tree), color = "red", fill = "black", shape = 23, size = 2, alpha = .8) +
  #geom_errorbar(aes(ymin=our_growthRate_VCF_mut_Tree_minusSD, ymax=our_growthRate_VCF_mut_Tree_plusSD), width=.2, color = "red")+
  geom_text(aes(label = paste0("n=", N_tips), y = 2.0))+
  ylab("Growth Rates (per year)") + #ylim(0, 2) + theme_bw() +
  ggtitle("Their growth rate (blue) vs. our growth rates using their trees (red square), \n ultrametric trees from VCF (red triangle), or mutation trees from VCF (red diamond with error bars +/- SD)") +
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


