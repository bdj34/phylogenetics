# Make PD6629 ultrametric from the non-ultrametric tree they give
library(ape)
library(phangorn)
library(ggplot2)
library(gridExtra)
library(ggtree)

rm(list = ls())

example_tree_PD6629 <- readRDS('~/Downloads/PD6629.RDS')

tree_with_mutations_theirs <- example_tree_PD6629[["pdx"]][["tree_ml"]]

agedf <- example_tree_PD6629$pdx$agedf
agedf$age <- agedf$age_at_sample_pcy
agedf <- agedf[,colnames(agedf) %in% c("tip.label", "age")]
tree_with_mutations_theirs$agedf <- agedf

tree_simple <- tree_with_mutations_theirs

############# Plot tree and check against the one in the paper
# Get sample and age of each tip
#sample_withTail <- sub(paste0(".*", pid), "", tree_simple$tip.label)
sample <- tree_simple$tip.label
sample[grep("ad", sample)] = "B"
sample[grep("d", sample)] = "A"
sample[-grep("A", sample)] = "B"
unique_samples <- c("A", "B")

n_samples <- length(unique_samples)
stopifnot(n_samples > 1 & n_samples < 4)

col_vec <- c()
age_vec <- c()
for ( l in 1:length(sample)) {
  if (sample[l] == unique_samples[1]) {
    col_vec[l] <- "red"
  } else if (sample[l] == unique_samples[2]) {
    col_vec[l] <- "blue"
  } else {
    col_vec[l] <- "green"
  }
}

sample[length(sample)] <- "zeros"
col_vec[length(sample)] <- "black"

# Add tip color info to the tree object
tree_simple$tip.color <- col_vec

# Plot (check against tree from paper) ###NOTE: One tip is colored wrong
print(ggtree(tree_simple, branch.length = tree_simple$edge.length) + 
        geom_tippoint(size = 2, color=tree_simple$tip.color)+
        layout_dendrogram())
plot_tree(tree_simple, cex.label = .5)


# Why is this not an integer to begin with? Adjusted for depth/sensitivity??
tree_simple$edge.length <- as.integer(round(tree_simple$edge.length))

# Make ultrametric
res=fit_tree(tree=tree_simple, model = "poisson_tree", switch_nodes = c(), niter = 1000)
ultraTree_from_theirMutTree <- res$ultratree

# Take their ultrametric, plot
all <- readRDS('~/Downloads/PDD_TELO.rds')
pid <- all[["PD6629"]]
patient_name <- pid$patient
tree <- pid$ultratree

#Plot ultra trees using plot_tree
plot_tree(tree, cex.label = 0)
plot_tree(res$ultratree, cex.label = 0)

# Plot ultra trees using ggtree (side by side)
theirs <- ggtree(tree)+layout_dendrogram()+ggtitle("Their ultrametric tree")
ours_using_theirs <- ggtree(ultraTree_from_theirMutTree)+layout_dendrogram()+ggtitle("Ultrametric tree from their mutation tree")
grid.arrange(theirs, ours_using_theirs, ncol = 2)


# Get our ultrametric tree from their VCF (see "mpboot_rtreefit_tree_construction_05102022.R")
pid <- "PD6629"
vcf_file <- paste0("~/rotation_fall2021/vcf_dir/", pid, '.vcf')

agedf <- all[[pid]]$agedf

# Read VCF and organize
vcf <- read.vcfR(vcf_file, verbose = FALSE )
data_mat_all <- vcf@gt[,2:dim(vcf@gt)[2]]
mut_info <- vcf@fix[,"INFO"]
row.names(data_mat_all) <- mut_info

ref_nuc <- vcf@fix[,"REF"]
alt_nuc <- vcf@fix[,"ALT"]

# Keep only SNVs for the tree building
rowsKeep <- which(nchar(ref_nuc) == 1 & nchar(alt_nuc) == 1)
rowsKeep_doubleCheck <- grep("SNV", row.names(data_mat_all))
rowsRemove <- grep("INDEL", row.names(data_mat_all))
stopifnot(all(rowsKeep == rowsKeep_doubleCheck))
stopifnot(length(rowsKeep) + length(rowsRemove) == dim(data_mat_all)[1])

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

# Read tree created with mpboot
filename <- paste0("~/rotation_fall2021/check_", pid, ".nexus.treefile")
# Load into R using ape
tree <- ape::read.tree(filename)

# Root the tree
tree <- root(tree, outgroup = "zeros", resolve.root = T)

# Assign arbitrary edge lengths
tree$edge.length <- c(0, rep(1, dim(tree$edge)[1]-1))
# Check that the tree looks normal
ggtree(tree, branch.length = tree$edge.length) + geom_tiplab(size = .4)

# Now we can assign all mutations to the trees (regardless of whether we called them originally)
mut <- matrix(data = NA, nrow = dim(data_mat_all)[1], ncol = dim(data_mat_all)[2]+1)
depth <- matrix(data = NA, nrow = dim(data_mat_all)[1], ncol = dim(data_mat_all)[2]+1)
colnames(mut) <- c("zeros", colnames(data_mat_all))
colnames(depth) <- c("zeros", colnames(data_mat_all))

for (j in 1:dim(data_mat_all)[1]){
  
  mut_entry <- as.numeric(sub(":.*", "", sub("^[^:]*:", '',data_mat_all[j,])))
  mut[j,] <- c(0, mut_entry)
  
  depth_entry <- mut_entry + as.numeric(sub(":.*", "", data_mat_all[j,]))
  depth[j,] <- c(10, depth_entry)
  
}

p.error <- c(0.000001, rep(0.01, dim(data_mat_all)[2]))
head(mut)
head(depth)

stopifnot(depth >= mut)

res_from_data = assign_to_tree(tree, mut, depth,error_rate=p.error)
tree_from_vcf <- res_from_data$tree
plot_tree(tree_from_vcf, cex.label = 1)

# Get age info
unique_ages <- unique(agedf$age_at_sample_pcy)
# Remove zeroes group
unique_ages_no_zeros <- sort(unique_ages[unique_ages > 1])

# Get sample and age of each tip
sample_withTail <- sub(paste0(".*", pid), "", tree_from_vcf$tip.label)
sample <- sub("_.*", "", sample_withTail)
sample[sample != "db"] = "B"
sample[sample == "db"] = "A"
unique_samples <- c("A", "B")

n_samples <- length(unique_samples)
stopifnot(n_samples > 1 & n_samples < 4)

col_vec <- c()
age_vec <- c()
for ( l in 1:length(sample)) {
  if (sample[l] == unique_samples[1]) {
    col_vec[l] <- "red"
    age_vec[l] <- unique_ages_no_zeros[1]
  } else if (sample[l] == unique_samples[2]) {
    col_vec[l] <- "blue"
    age_vec[l] <- unique_ages_no_zeros[2]
  } else {
    col_vec[l] <- "green"
    age_vec[l] <- unique_ages_no_zeros[3]
  }
}

sample[1] <- "zeros"
col_vec[1] <- "black"
age_vec[1] <- unique_ages[unique_ages < 1]

# Add tip color info to the tree object
tree_from_vcf$tip.color <- col_vec

# Now we want the number of mutations from root (zeros) to each tip
dist_nodes_mat <- dist.nodes(tree_from_vcf)

# Convert to ultrametric
age_vec_to_df <- data.frame("tip.label" = tree_from_vcf$tip.label, "age" = age_vec)
tree_from_vcf$agedf <- age_vec_to_df

res_ultra_from_VCF=fit_tree(tree=tree_from_vcf, model = "poisson_tree", switch_nodes = c(), niter = 1000)
our_tree_vcf <- res_ultra_from_VCF$ultratree

our_tree_plot_from_vcf <- ggtree(our_tree_vcf)+layout_dendrogram()+ggtitle("Ultrametric Tree From VCF")

grid.arrange(theirs, ours_using_theirs, our_tree_plot_from_vcf, ncol = 3)

