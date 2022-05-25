# Get age/rate estimates for all samples
library(ape)
library(phangorn)
library(treemut)
library(rtreefit)
library(ggtree)
library(ggplot2)
library(vcfR)
library(castor)
library(gridExtra)

rm(list = ls())

# Set threshold number of tips at which we will calculate clade age (default = 10)
threshold_n_tips <- 10

all <- readRDS('~/Downloads/PDD_TELO.rds')
clade_info <- read.table("~/rotation_fall2021/clade_info.txt", header = T)
clade_info$lambda <- clade_info$lambda/100
clade_info$lower_lambda <- clade_info$lower_lambda/100
clade_info$upper_lambda <- clade_info$upper_lambda/100

# Make df for saving
df <- clade_info
df$our_clade_age <- NA
df$their_clade_age <- NA
df$fit_clade_lambda <- NA
df$N_tips <- NA

for (pid_num in 1:length(all)) {
  #pid_num <- 2
  
  pid <- all[[pid_num]]
  
  patient_name <- pid$patient
  
  tree <- pid$ultratree
  meta <- pid$agedf
  meta$driver3 <- gsub(",", ":", meta$driver3)
  
  # Ignore 1q+ in PD5182 (not in table from paper)
  if(patient_name == "PD5182") {
    meta$driver3 <- gsub("1q\\+:9pUPD:JAK2", "JAK2:9pUPD", meta$driver3)
  }
  
  # Go through each driver at each timepoint
  for (driver in unique(meta$driver3)) {
    for (age in unique(meta$age_at_sample_pcy)) {
      
      count <- sum(meta$driver3 == driver & meta$age_at_sample_pcy == age)
      
      if (count > threshold_n_tips & driver != "WT") {
        #print("Threshold exceeded")
        
        # Get indexes of the samples that have this age and driver combo
        rowsKeep <- which(meta$driver3 == driver & meta$age_at_sample_pcy == age)
        tipsKeep <- meta$tip.label[rowsKeep]
        
        subset_tree <- get_subtree_with_tips(tree, only_tips = tipsKeep, force_keep_root = T)
        
        print(ggtree(subset_tree$subtree) + ggtitle(paste0(round(age,2), " yrs with ", driver)) + layout_dendrogram())
        
        # Have subset tree, now get internal and external branch lengths
        ## first get the node numbers of the tips
        nodes<-sapply(subset_tree$subtree$tip.label,function(x,y) which(y==x),y=subset_tree$subtree$tip.label)
        ## then get the edge lengths for those nodes
        ext_lengths<-setNames(subset_tree$subtree$edge.length[sapply(nodes,
                                                                         function(x,y) which(y==x),y=subset_tree$subtree$edge[,2])],names(nodes))
        
        ggtree(subset_tree$subtree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
        
        # Find the node that only has one out node (root) and take that out node as the new root (closest to where mut occurs)
        root <- as.numeric(names(which(table(subset_tree$subtree$edge[,1])==1))[1])
        new_root <- subset_tree$subtree$edge[which(subset_tree$subtree$edge[,1]==root),2]
        total <- dist.nodes(subset_tree$subtree)[1:length(subset_tree$subtree$tip.label),new_root]
        
        # Because we have ultrametric trees, total edge length from root to tip is same for all tips
        #stopifnot(total[1] == total[2] & total[5] == total[9])
        
        
        internal <- total-ext_lengths
        names(internal) <- subset_tree$subtree$tip.label
        
        # Calculate age. Note: lengths ar already in units of mutations
        n <- length(subset_tree$subtree$tip.label)
        
        # Determine if driver is in the inputted table
        bool <- F
        for (l in 1:length(clade_info$driver)) {
          bool <- all(unlist(strsplit(driver, ":")) %in% unlist(strsplit(clade_info$driver[l], ":"))) &
            clade_info$pid[l] == patient_name
          if (bool) {
            row <- l
            break
          }            
        }
        
        
        if (bool) {
          #print("driver in clade_info$driver")
          lambda <- log(1+clade_info$lambda[row])
          clade_age_postConc <- clade_info$clade_age_PC[row]
          
          ourAge <- ((n-1)*sum(ext_lengths))/(n*lambda*sum(internal)) + (1/lambda)*(1+sum(1/c(1:(n-1))))
          theirAge <- age - clade_age_postConc
          our_lambda_their_age <- (((n-1)*sum(ext_lengths))/(n*sum(internal)) + (1+sum(1/c(1:(n-1)))))/theirAge
          
          df$our_clade_age[row] <- ourAge
          df$their_clade_age[row] <- theirAge
          df$fit_clade_lambda[row] <- our_lambda_their_age
          df$N_tips[row] <- count
          
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
write.csv(out.df, "~/rotation_fall2021/clade_estimates_04202022.csv")

ggplot(out.df, aes(x = pid)) + geom_point(aes(y = their_clade_age), color = "blue") + 
  geom_point(aes(y = our_clade_age), color = "red") + ylab("Clade age (yr)")

       
       
