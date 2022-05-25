# Age and site freq from ultrametric trees (double checking)
library(ape)
library(phangorn)
library(treemut)
library(rtreefit)
library(ggtree)
library(ggplot2)
library(vcfR)
library(castor)
library(gridExtra)
library(adephylo)

rm(list = ls())

# Set threshold number of tips at which we will calculate clade age (default = 10)
threshold_n_tips <- 10

all <- readRDS('~/Downloads/PDD_TELO.rds')
clade_info <- read.table("~/rotation_fall2021/clade_info.txt", header = T)
clade_info$lambda <- log(clade_info$lambda/100 + 1)
clade_info$lower_lambda <- log(clade_info$lower_lambda/100 + 1)
clade_info$upper_lambda <- log(clade_info$upper_lambda/100 + 1)

# Make df for saving
df <- data.frame("pid" = NA, "driver" = NA, "lambda" = NA, "lower_lambda" = NA, 
                 "upper_lambda" = NA, "clade_age_PC" = NA, "lower_clade_age_PC" = NA, 
                 "upper_clade_age_PC" = NA, "fit_clade_lambda" = NA, "our_clade_age" = NA, 
                 "clade_age_ABC" = NA, "clade_age_lower_ABC" = NA, "clade_age_upper_ABC" = NA,
                 "clade_age_tree" = NA,  "N_tips" = NA, 
                 "internal_edge_length" = NA, "external_edge_length" = NA, 
                 "age_at_sample" = NA)

count <- 1
for (pid_num in 1:length(all)) {
  
  #pid_num <- 2
  
  pid <- all[[pid_num]]
  
  patient_name <- pid$patient
  
  tree <- pid$ultratree
  meta <- pid$agedf
  meta$driver3 <- gsub(",", ":", meta$driver3)
  
  # Ignore 1q+ in PD5182 (not in table from paper) and correct first sample by adding 9pUPD
  # Note: PD5182 is weird and poorly annotated in the given df so ignore if necessary
  if(patient_name == "PD5182") {
    meta$driver3 <- gsub("1q\\+:9pUPD:JAK2", "JAK2:9pUPD", meta$driver3)
    meta$driver3 <- gsub("^JAK2$", "JAK2:9pUPD", meta$driver3)
    meta$driver3[meta$tip.label == "db19"] <- "JAK2"
    meta$driver3[meta$driver3 == "JAK2:9pUPD" & meta$age_at_sample_pcy > 40] <- "JAK2:9pUPD:1q+"
    meta$driver3[meta$tip.label == "dl05"] <- "JAK2:9pUPD"
    next # Just skip for now
  }
  
  
  
  # Go through each driver at each timepoint
  for (driver in unique(meta$driver3)) {
    for (age in unique(meta$age_at_sample_pcy)) {
      #driver <- unique(meta$driver3)[2]
      #age <- unique(meta$age_at_sample_pcy)[1]
      
      n <- sum(meta$driver3 == driver & meta$age_at_sample_pcy == age)
      
      if (n >= threshold_n_tips & driver != "WT") {
        
        # Get indexes of the samples that have this age and driver combo
        rowsKeep <- which(meta$driver3 == driver & meta$age_at_sample_pcy == age)
        tipsKeep <- meta$tip.label[rowsKeep]
        
        subset_tree <- get_subtree_with_tips(tree, only_tips = tipsKeep, force_keep_root = T)
        subtree <- subset_tree$subtree
        
        #print(ggtree(subtree) + ggtitle(paste0(round(age,2), " yrs with ", driver)) + layout_dendrogram())
        
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
        
        write.csv(descendant_df, paste0("~/rotation_fall2021/site_freq_plots/descendant_info_", 
                                       patient_name, "_age_", round(age,2), "_driver_", gsub(":", "_", driver), ".csv"))
        
        # Finally, sum up and record the edge lengths which are shared by the same number of cells
        site_freq <- data.frame("n_cells" = unique(descendant_df$n_cells), "freq" = NA)
        for (k in site_freq$n_cells) {
          site_freq$freq[site_freq$n_cells == k] <- sum(descendant_df$Edge_length[descendant_df$n_cells == k])
        }
        
        write.csv(site_freq, paste0("~/rotation_fall2021/site_freq_plots/site_freq_", 
                                       patient_name, "_age_", round(age,2), "_driver_", gsub(":", "_", driver), ".csv"))
        
        # Write a function to plot alongside data for site freq spectrum which fits 1/k^2 based on frequency of muts in 2 cells
        constant <- 4*site_freq$freq[site_freq$n_cells == 2]
        site_freq_fit <- function(x) {
          return(constant/(x^2))
        }
        
        # Optionally plot and save
        #pdf(paste0("~/rotation_fall2021/site_freq_plots/", patient_name, "_age_", round(age,2), "_driver_", gsub(":", "_", driver), ".pdf"))
        #print(ggplot(site_freq) + geom_point(aes(x = n_cells, y = freq)) + geom_function(fun = site_freq_fit) +
        #        ggtitle(paste0("Site frequency spectrum ", patient_name, " sample age: ", round(age,2), " driver: ", driver)))
        #dev.off()
        
        # The sum of edge lengths in site_freq dataframe should equal the total internal branch lengths
        site_freq_internal <- sum(site_freq$freq)
        
        # From subtree, get  external branch lengths
        ## first get the node numbers of the tips
        tips <- sapply(subtree$tip.label, function(x,y) which(y==x), y = subtree$tip.label)
        ## then get the edge lengths for edges ending at the tips
        ext_lengths <- setNames(subtree$edge.length[sapply(tips, function(x,y) which(y==x), y = subtree$edge[,2])], names(tips))
        
        #ggtree(subtree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
        #plot_tree(subtree)
        
        ####### Calculate age. Note: lengths are in time units (ratio; shouldn't matter)
        # First get n, number of sampled cells in this lineage
        n <- length(subtree$tip.label)
        
        # Distance from root to tip gives us their tree age estimate
        new_root <- length(subtree$tip.label) + 2
        clade_age_tree <- dist.nodes(subset_tree$subtree)[1:length(subset_tree$subtree$tip.label),new_root][1]
        
        # Determine if driver is in the inputted clade_info table (need their growth rate)
        bool <- F
        for (l in 1:length(clade_info$driver)) {
          bool <- all(unlist(strsplit(driver, ":")) %in% unlist(strsplit(clade_info$driver[l], ":"))) &
            all(unlist(strsplit(clade_info$driver[l], ":")) %in% unlist(strsplit(driver, ":"))) &
            clade_info$pid[l] == patient_name
          if (bool) {
            row <- l
            break
          }            
        }
        
        if (bool) {
          #print("driver in clade_info$driver")
          lambda <- clade_info$lambda[row]
          clade_age_ABC_PC <- clade_info$clade_age_PC[row]
          clade_age_lower_ABC_PC <- clade_info$lower_clade_age_PC[row]
          clade_age__upper_ABC_PC <- clade_info$upper_clade_age_PC[row]
          
          ourAge <- ((n-1)*sum(ext_lengths))/(n*lambda*site_freq_internal) + (1/lambda)*(1+sum(1/c(1:(n-1))))
          theirAge_ABC <- age - clade_age_ABC_PC
          theirAge_upper_ABC <- age - clade_age_lower_ABC_PC
          theirAge_lower_ABC <- age - clade_age__upper_ABC_PC
          our_lambda_their_age <- (((n-1)*sum(ext_lengths))/(n*site_freq_internal) + (1+sum(1/c(1:(n-1)))))/clade_age_tree
          
          df[count,] <- c(clade_info[row,], our_lambda_their_age, ourAge, theirAge_ABC, 
                          theirAge_lower_ABC, theirAge_upper_ABC,
                          clade_age_tree, n, site_freq_internal, sum(ext_lengths), age)
          
          count <- count + 1
          #print(paste0("Our tumor age: ", ourAge, "   their tumor age: ", theirAge, "   their tumor lambda: ", lambda, "   fit lambda: ", our_lambda_their_age))
          
        }else {
          pid_clades <- clade_info[clade_info$pid == patient_name,]
          print(paste0("In patient ", patient_name, ", ", driver, " not in ", pid_clades$driver))
        }
      }
    }
  }
}

out.df <- df[!is.na(df$fit_clade_lambda),]
write.csv(out.df, "~/rotation_fall2021/clade_estimates_04282022.csv")

out.df$driver_and_pid_and_age <- paste0(out.df$pid, " ", out.df$driver, " ", round(out.df$age_at_sample, 2))
out.df$driver_and_pid <- paste0(out.df$pid, " ", out.df$driver)

ggplot(out.df, aes(x = driver_and_pid_and_age)) + geom_point(aes(y = clade_age_ABC), color = "blue") + 
  geom_errorbar(aes(ymin=clade_age_lower_ABC, ymax=clade_age_upper_ABC), width=.2, color = "blue") +
  geom_point(aes(y = our_clade_age), color = "red") + 
  geom_point(aes(y = clade_age_tree), color = "green") + ylab("Clade age (yr)") + 
  ggtitle("Their ABC age (blue), their tree age (green) and our age using their growth rate (red)") +
  theme(axis.text.x = element_text(face = "bold", angle = 70, vjust = 1, hjust=1)) 

ggplot(out.df, aes(x = driver_and_pid)) + geom_point(aes(y = lambda), color = "blue") + 
  geom_errorbar(aes(ymin=lower_lambda, ymax=upper_lambda), width=.2, color = "blue") + 
  geom_point(aes(y = fit_clade_lambda), color = "red") + ylab("Lambda") + ylim(0, 2) +
  ggtitle("Their lambda (blue) vs. what our lambda would be to match tree age (red)") +
  theme(axis.text.x = element_text(face = "bold", angle = 70, vjust = 1, hjust=1))
