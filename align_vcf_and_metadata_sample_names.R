# Align VCF sample names with agedf sample names
library(vcfR)
library(readr)

rm(list = ls())

outDir <- "~/rotation_fall2021/vcfs_with_matched_sample_names/"
dir.create(outDir)

# Import information about samples including drivers, ultrametric trees from paper, etc.
all_ultra <- readRDS('~/Downloads/PDD_TELO.rds')

pids <- names(all_ultra)

for (i in 1:length(pids)) {

  patientName <- pids[i]
  patientInfo <- all_ultra[[i]]
  
  vcf_file <- paste0("~/rotation_fall2021/vcf_dir/", patientName, ".vcf")
  
  # Read VCF and organize
  vcf <- read.vcfR(vcf_file, verbose = FALSE )
  data_mat_all <- vcf@gt[,2:dim(vcf@gt)[2]]
  
  # Check that the first column of the vcf is format
  stopifnot(colnames(vcf@gt)[1] == "FORMAT")
  
  sampleNamesLong_vcf <- colnames(data_mat_all)
  sampleNamesNoPID_vcf <- gsub(patientName, "", sampleNamesLong_vcf)
  
  # Get agedf (metadata) and make sure zeros is last tip
  agedf <- patientInfo$agedf
  stopifnot(agedf$tip.label[length(agedf$tip.label)] == "zeros")
  samples_from_agedf <- agedf$tip.label
  n_samples <- length(unique(agedf$age_at_sample_pcy[agedf$age_at_sample_pcy > 1]))
  
  vcf_correctName <- vcf
  
  # Match agedf metadata to vcf data by parsing numbers and matching (only works for single sample patients rn)
  if (n_samples == 1) {
    
    if (patientName %in% c("PD5179", "PD9478", "PD5117")) {
      
      agedf_col_ordered_by_vcf <- match(sampleNamesNoPID_vcf, samples_from_agedf)
      
    }else {
      
      sampleNumber_vcf <- parse_number(sampleNamesNoPID_vcf)
      sampleNumber_agedf <- parse_number(agedf$tip.label[1:(length(agedf$tip.label)-1)])
      
      agedf_col_ordered_by_vcf <- match(sampleNumber_vcf, sampleNumber_agedf)
    }
    
    stopifnot(all(!is.na(agedf_col_ordered_by_vcf)))
    colnames(vcf_correctName@gt) <- c("FORMAT", samples_from_agedf[agedf_col_ordered_by_vcf])
    
    
    
  } else {
    if (patientName %in% c("PD5182", "PD4781")) {
      samples_vcf <- gsub("_lo00", "", sampleNamesNoPID_vcf)
      agedf_col_ordered_by_vcf <- match(samples_vcf, samples_from_agedf)
    } else if (patientName == "PD6629") {
      samples_vcf <- gsub("b_lo00", "", sampleNamesNoPID_vcf)
      agedf_col_ordered_by_vcf <- match(samples_vcf, samples_from_agedf)
    } else if (patientName == "PD6646") {
      samples_vcf <- gsub("d_lo0", "", sampleNamesNoPID_vcf)
      samples_vcf <- gsub("m_lo00", "", samples_vcf)
      samples_vcf[grep("m", sampleNamesNoPID_vcf)] <- paste0(samples_vcf[grep("m", sampleNamesNoPID_vcf)], "m")
      samples_vcf[grep("d", sampleNamesNoPID_vcf)] <- paste0(samples_vcf[grep("d", sampleNamesNoPID_vcf)], "d")
      agedf_col_ordered_by_vcf <- match(samples_vcf, samples_from_agedf)
    }
    
    
    stopifnot(all(!is.na(agedf_col_ordered_by_vcf)))
    colnames(vcf_correctName@gt) <- c("FORMAT", samples_from_agedf[agedf_col_ordered_by_vcf])
  }
  
  saveRDS(vcf_correctName, paste0(outDir, patientName, ".rds"))

  
}


