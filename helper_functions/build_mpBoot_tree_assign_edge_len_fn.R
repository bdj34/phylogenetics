build_mpBoot_tree_assign_edge_len_fn <- function(vcf, dir, patientName, use_saved_tree = T) {
  
  outFile <- paste0(dir, patientName, "_tree_with_mutations_from_VCF.rds")
  if (file.exists(outFile) & use_saved_tree) {
    tree <- readRDS(outFile)
    return (tree)
  }
  
  data_mat_all <- vcf@gt[,2:dim(vcf@gt)[2]]
  mut_info <- vcf@fix[,"INFO"]
  row.names(data_mat_all) <- mut_info
  
  ref_nuc <- vcf@fix[,"REF"]
  alt_nuc <- vcf@fix[,"ALT"]
  
  # Excluded region = CNA regions
  is_in_excluded_region <- as.numeric(gsub(".*REGION=", "", vcf@fix[,"INFO"]))
  
  # Keep only SNVs OUTSIDE EXCLUDED REGION! for the tree building
  rowsKeep <- which(nchar(ref_nuc) == 1 & nchar(alt_nuc) == 1 & is_in_excluded_region == 0)
  
  data_mat_cut <- data_mat_all[rowsKeep,]
  
  ref_nuc_cut <- ref_nuc[rowsKeep]
  alt_nuc_cut <- alt_nuc[rowsKeep]
  
  #Make numeric matrix with 0 (not mutated) or 1 (mutated)
  data_matrix_binarize <- matrix(0, nrow = dim(data_mat_cut)[1], ncol = dim(data_mat_cut)[2])
  for (j in 1:dim(data_mat_cut)[1]){
    data_matrix_binarize[j,] <- suppressWarnings(as.numeric(sub('^([^:]+:[^:]+):', '', data_mat_cut[j,])))
  }
  row.names(data_matrix_binarize) <- mut_info[rowsKeep]
  colnames(data_matrix_binarize) <- colnames(data_mat_all)
  
  # Substitute unknown variants with the correct IUPAC system name
  return_unknown_nuc_symbol <- function(ref, alt) {
    vec <- c(ref, alt)
    if (all(vec %in% c("A", "T"))) {
      return("W")
    } else if (all(vec %in% c("C", "G"))){
      return("S")
    } else if (all(vec %in% c("A", "C"))){
      return("M")
    } else if (all(vec %in% c("G", "T"))){
      return("K")
    } else if (all(vec %in% c("A", "G"))){
      return("R")
    } else if (all(vec %in% c("C", "T"))){
      return("Y")
    } else {
      stop("error")
    }
  }
  
  nexus_matrix <- matrix("error", nrow = dim(data_matrix_binarize)[1], ncol = dim(data_matrix_binarize)[2]+1)
  for (j in 1:dim(nexus_matrix)[1]) {
    zero_cols <- which(data_matrix_binarize [j,] == 0)
    one_cols <- which(data_matrix_binarize [j,] == 1)
    unknown_cols <- which(is.na(data_matrix_binarize[j,]))
    
    nexus_matrix[j,one_cols] <- alt_nuc_cut[j]
    nexus_matrix[j,zero_cols] <- ref_nuc_cut[j]
    nexus_matrix[j,unknown_cols] <- return_unknown_nuc_symbol(ref_nuc_cut[j], alt_nuc_cut[j])
  }
  
  colnames(nexus_matrix) <- c(colnames(data_mat_all), "zeros")
  nexus_matrix[,"zeros"] <- ref_nuc_cut
  
  # Write nexus file out to make tree
  file <- paste0(dir, '/nexus_matrix_mpboot_input_', patientName, '.nexus')
  write.nexus.data(t(nexus_matrix), file, format = "dna", datablock = TRUE,
                   interleaved = F, charsperline = NULL,
                   gap = NULL, missing = NULL)
  
  # Run mpboot using system
  system(paste0("~/mpboot-avx-1.1.0-MacOSX/bin/mpboot -s ", file, " -st DNA"))
  
  # Read tree created with mpboot
  filename <- paste0(tree_dir, "/nexus_matrix_mpboot_input_",
                     patientName, ".nexus.treefile")
  # Load into R using ape
  tree <- ape::read.tree(filename)
  
  # Root the tree
  tree <- root(tree, outgroup = "zeros", resolve.root = T)
  
  # Assign arbitrary edge lengths
  tree$edge.length <- c(rep(1, dim(tree$edge)[1]-1), 0)
  # Check that the tree looks normal
  ggtree(tree, branch.length = tree$edge.length) + geom_tiplab(size = 3)
  
  # Now we can assign all mutations to the trees (regardless of whether we called them originally)
  mut <- matrix(data = NA, nrow = dim(data_mat_all)[1], ncol = dim(data_mat_all)[2]+1)
  depth <- matrix(data = NA, nrow = dim(data_mat_all)[1], ncol = dim(data_mat_all)[2]+1)
  colnames(mut) <- c(colnames(data_mat_all), "zeros")
  colnames(depth) <- c(colnames(data_mat_all), "zeros")
  
  for (j in 1:dim(data_mat_all)[1]){
    
    mut_entry <- as.numeric(sub(":.*", "", sub("^[^:]*:", '',data_mat_all[j,])))
    mut[j,] <- c(mut_entry, 0)
    
    depth_entry <- mut_entry + as.numeric(sub(":.*", "", data_mat_all[j,]))
    depth[j,] <- c(depth_entry, 10)
    
  }
  
  # Sample specific????
  p.error <- c(rep(0.01, dim(data_mat_all)[2]), 0.000001)
  #head(mut)
  #head(depth)
  
  stopifnot(depth >= mut)
  
  res_from_data = assign_to_tree(tree, mut, depth,error_rate=p.error)
  tree_with_mutations <- res_from_data$tree
  
  saveRDS(tree_with_mutations, outFile)
  
  return(tree_with_mutations)
}
