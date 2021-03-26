#!/usr/bin/env Rscript

cl <- commandArgs(trailingOnly = TRUE)

marker_file_1 <- cl[1]
species_1 <- cl[2]
marker_file_2 <- cl[3]
species_2 <- cl[4]
ortholog_mapping_file <- cl[5]
pval_limit <- as.numeric(cl[6])
ngenes <- as.numeric(cl[7])
outfile <- cl[8]

min_overlap <- 0

if (ngenes == 'all'){
    ngenes <- 1000000
}

markers_1 <- subset(read.delim(marker_file_1, stringsAsFactors = FALSE), pvals_adj < pval_limit & rank < ngenes)
markers_2 <- subset(read.delim(marker_file_2, stringsAsFactors = FALSE), pvals_adj < pval_limit & rank < ngenes)

ortholog_mapping <- read.delim(ortholog_mapping_file)
species_1_orth_col <- paste0(species_1, '_gene_id')
species_2_orth_col <- paste0(species_2, '_gene_id')

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

if (length(valid_intersects) == 0){
  write(paste0("No cell groups overlap between ", marker_file_1, ' and ', marker_file_2, ' at a proportion of at least ', min_overlap), stderr())
  q(status = 1)
}

cluster_matches <- data.frame(do.call(rbind, unlist(lapply(names(valid_intersects), function(x) strsplit(x, '\\.')), recursive = FALSE)))
colnames(cluster_matches) <- c(paste(species_1, 'cluster'), paste(species_2, 'cluster'))
cluster_matches$intersect_prop <- unlist(gene_intersect_props[names(valid_intersects)])
cluster_matches$intersect <- unlist(lapply(gene_intersects[names(valid_intersects)], function(x) length(x)))
cluster_matches$intersect_gene_ids <- unlist(lapply(gene_intersects[names(valid_intersects)], function(x) paste(x, collapse = ', ')))
cluster_matches$intersect_gene_symbols <- unlist(lapply(gene_intersects[names(valid_intersects)], function(x) paste(ortholog_mapping[match(x, ortholog_mapping[[species_2_orth_col]]),paste0(species_2, '_gene_name')], collapse = ', ')))

# Get the set of matching clusters passing criteria for each query
cluster_matches <- cluster_matches[order(cluster_matches$intersect_prop, decreasing = TRUE),]
cluster_matches_by_query <- split(cluster_matches, cluster_matches[[paste(species_1, 'cluster')]])

# Order such that strongest matches are first 
match_order <- order(unlist(lapply(cluster_matches_by_query, function(x) max(x$intersect_prop))), decreasing = TRUE)

# Re-bind the results for output

write.table(do.call(rbind, cluster_matches_by_query[match_order]), file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
