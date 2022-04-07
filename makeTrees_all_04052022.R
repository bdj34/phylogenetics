# Make all trees and color by sample.
# Color of nodes has been manually checked by Brian and should be good.
# Note: This uses their variant calling and excludes "?"'s from the tree building
library(vcfR)
library(ggplot2)
rm(list =ls())
setwd('~/rotation_fall2021/vcf_files')

vcf_files <- list.files(pattern = '.vcf')

for (i in 1:length(vcf_files)) {

  # Get the sample name
  pid <- sub(".vcf", "", vcf_files[i])
  
  # Read VCF and organize
  vcf <- read.vcfR( vcf_files[i], verbose = FALSE )
  data_mat_all <- vcf@gt[,2:dim(vcf@gt)[2]]
  mut_info <- vcf@fix[,"INFO"]
  row.names(data_mat_all) <- mut_info
  
  #Make numeric matrix with 0 (not mutated) or 1 (mutated)
  data_matrix_binarize <- matrix(0, nrow = dim(data_mat_all)[1], ncol = dim(data_mat_all)[2])
  for (j in 1:dim(data_mat_all)[1]){
    data_matrix_binarize[j,] <- suppressWarnings(as.numeric(sub('^([^:]+:[^:]+):', '', data_mat_all[j,])))
  }
  row.names(data_matrix_binarize) <- mut_info
  colnames(data_matrix_binarize) <- colnames(data_mat_all)
  
  # Remove variants that have NA's (probably a better way to do this)
  data_matrix_noNA_all <- data_matrix_binarize[!is.na(rowSums(data_matrix_binarize)),]
  
  # Remove variants that are mutated in all cells or no cells (these may already be removed?)
  data_matrix_noNA <- data_matrix_noNA_all[rowSums(data_matrix_noNA_all) != 0 &
                                             rowSums(data_matrix_noNA_all) != dim(data_matrix_noNA_all)[2], ]
  
  # Convert our matrix to a phangorn object
  phyDat <- MatrixToPhyDat(t(data_matrix_noNA))
  
  # Use phangorn parsimony ratchet function to make the tree (with no edge lengths)
  treeRatchet <- phangorn::pratchet(phyDat,
                                    maxit = 10^6, minit = 100, k = 100, method = "fitch",
                                    rearrangements = "SPR", trace = 1)
  
  ggtree(treeRatchet) + layout_dendrogram()
  
  # Set edge lengths using acctran function
  treeRatchet <- acctran(treeRatchet, phyDat)
  treeRatchet_root <- midpoint(treeRatchet)
  
  #ggtree(treeRatchet_root) + layout_dendrogram()
  
  #sort(treeRatchet_root$tip.label)
  
  # Get sample of each tip
  sample_withTail <- sub(paste0(".*", pid), "", treeRatchet_root$tip.label)
  sample <- sub("_.*", "", sample_withTail)
  if (pid == "PD6629") {
    sample[sample != "db"] = "A"
    sample[sample == "db"] = "B"
    unique_samples <- c("B", "A")
  }else if (pid == "PD5182"){
    sample[sample != "db" & sample != "dl"] = "A"
    sample[sample == "db"] = "B"
    sample[sample == "dl"] = "C"
    unique_samples <- c("B", "A", "C")
  } else if (pid == "PD6646") {
    unique_samples <- c("m", "d")
  }else {
    unique_samples <- unique(sample)
  }

  n_samples <- length(unique_samples)
  if (all(sample == sample_withTail)) {
    n_samples <- 1
    unique_samples <- c("A")
    sample <- rep("A", length(sample))
  }
  
  col_vec <- c()
  for ( l in 1:length(sample)) {
    if (sample[l] == unique_samples[1]) {
      col_vec[l] <- "red"
    } else if (sample[l] == unique_samples[2]) {
      col_vec[l] <- "blue"
    } else {
      col_vec[l] <- "green"
    }
  }
  
  # Add tip color info to the tree object
  treeRatchet_root$tip.color <- col_vec
  
  p <- ggtree(treeRatchet_root) + geom_tippoint(size = 2, color=treeRatchet_root$tip.color)+
    layout_dendrogram()
  
  pdf(paste0('~/rotation_fall2021/makeTrees_04052022/', "patient_", pid, '.pdf'))
  print(p)
  dev.off()
  
  saveRDS(treeRatchet_root, paste0('~/rotation_fall2021/makeTrees_04052022/', "patient_", pid, '.rds'))

}

