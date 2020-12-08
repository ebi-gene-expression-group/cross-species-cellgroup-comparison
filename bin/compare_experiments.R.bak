cl <- commandArgs(trailingOnly = TRUE)

experiment_1 <- cl[1]
species_1 <- cl[2]
experiment_2 <- cl[3]
species_2 <- cl[4]
ortholog_mapping_dir <- cl[5]
scxa_results_dir <- cl[6]

# Read ortholog mapping file and check for species

ortholog_mapping_file = file.path(ortholog_mapping_dir, paste(paste(species_1, species_2, sep='_'), 'txt', sep='.'))

if (! file.exists(ortholog_mapping_file)){
  stop(paste('No ortholog mapping file at', ortholog_mapping_file))
}

ortholog_mapping <- read.delim(ortholog_mapping_file)
species_1_orth_col <- paste0(species_1, '_gene_id')
species_2_orth_col <- paste0(species_2, '_gene_id')

# Find marker files in both experiments

bundle_dir_1 <- file.path(scxa_results_dir, experiment_1, species_1, 'bundle')
bundle_dir_2 <- file.path(scxa_results_dir, experiment_2, species_2, 'bundle')

manifest_file_1 <- file.path(scxa_results_dir, experiment_1, species_1, 'bundle', 'MANIFEST')
manifest_file_2 <- file.path(scxa_results_dir, experiment_2, species_2, 'bundle', 'MANIFEST')

manifest_1 <- read.delim(manifest_file_1, sep="\t", header = FALSE, stringsAsFactors = FALSE)
markers_manifest_1 <- subset(manifest_1, V1 %in% c('cluster_markers', 'meta_markers'))

manifest_2 <- read.delim(manifest_file_2, sep="\t", header = FALSE, stringsAsFactors = FALSE)
markers_manifest_2 <- subset(manifest_2, V1 %in% c('cluster_markers', 'meta_markers'))

marker_files_1 <- structure(file.path(bundle_dir_1, markers_manifest_1$V2), names = markers_manifest_1$V3)
marker_files_2 <- structure(file.path(bundle_dir_2, markers_manifest_2$V2), names = markers_manifest_2$V3)

# Function to compare two marker tables via ortholog identifiers

compare_marker_sets <- function(marker_file_1, marker_file_2, ortholog_mapping, species_1, species_2, species_1_orth_col, species_2_orth_col, pval_limit = 0.05, min_overlap = 0.1){
  markers_1 <- subset(read.delim(marker_file_1, stringsAsFactors = FALSE), pvals_adj < pval_limit)
  markers_2 <- subset(read.delim(marker_file_2, stringsAsFactors = FALSE), pvals_adj < pval_limit)
  
  # Sort by p value 
  
  markers_1 <- markers_1[order(markers_1$pvals_adj), ]
  markers_2 <- markers_2[order(markers_2$pvals_adj), ]
  
  # Add a column mapping the genes of the first dataset to the genes of the second
  markers_1$orth_id <- ortholog_mapping[match(markers_1$genes, ortholog_mapping[[species_1_orth_col]]), species_2_orth_col]
  
  # Split the markers tables by cluster and compare across species
  cluster_markers_1 <- split(markers_1, markers_1$cluster)
  cluster_markers_2 <- split(markers_2, markers_2$cluster)
  
  gene_intersects <- unlist(lapply(cluster_markers_1, function(cm1){
    lapply(cluster_markers_2, function(cm2){
      intersect(cm1$orth_id, cm2$genes) 
    })
  }), recursive = FALSE)
  
  gene_intersect_props <- unlist(lapply(cluster_markers_1, function(cm1){
    lapply(cluster_markers_2, function(cm2){
      length(intersect(cm1$orth_id, cm2$genes)) / length(which(! is.na(cm1$orth_id)))
    })
  }), recursive = FALSE)
  
  valid_intersects <- which(unlist(lapply(gene_intersect_props, function(x) x > min_overlap)))
  
  cluster_matches <- data.frame(do.call(rbind, unlist(lapply(names(valid_intersects), function(x) strsplit(x, '\\.')), recursive = FALSE)))
  colnames(cluster_matches) <- c(paste(species_1, 'cluster'), paste(species_2, 'cluster'))
  cluster_matches$intersect <- unlist(lapply(gene_intersects[names(valid_intersects)], function(x) length(x)))
  cluster_matches$intersect_prop <- unlist(gene_intersect_props[names(valid_intersects)])

  list(cluster_matches, gene_intersects)
}

marker_files_1 <- list.files("E-MTAB-5061_markers/", full.names = TRUE)
marker_files_2 <- list.files("E-MTAB-8483_markers/", full.names = TRUE)

lapply(marker_files_2, function(x){
  foo <- compare_marker_sets(marker_files_1[1], x, ortholog_mapping = ortholog_mapping, species_1 = species_1, species_2 = species_2, species_1_orth_col = species_1_orth_col, species_2_orth_col = species_2_orth_col)
  foo[[1]]
})


print(paste('Marker files 1:'))
print(marker_files_1)
print(paste('Marker files 2:'))
print(marker_files_2)
#print(paste('Marker files 2:', paste(marker_files_2, collapse = ',')))