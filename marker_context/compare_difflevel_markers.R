#!/usr/bin/env Rscript

# Take two marker gene lists derived from different contexts, and see how the
# overlap in genes changes as we make p value more and more stringent.

cl <- commmandArgs(trailingOnly = TRUE)

marker_file_1 <- "~/output/E-ENAD-15.pancreas.markers.tsv"
marker_file_2 <- "~/output/E-ENAD-15.markers.tsv"

markers1 <- read.delim(marker_file_1, stringsAsFactors = FALSE)
markers2 <- read.delim(marker_file_2, stringsAsFactors = FALSE)

common_celltypes <- unique(intersect(markers1$cluster, markers2$cluster))

# Find the proporition of overlapping marker genes at progressivlely relaxed
# thrsholds, starting at the minimum adjusted p value for each cell type

prop_overlaps <- do.call(rbind, lapply(common_celltypes, function(ct){
  min_pval <- floor(log10(min(markers1$pvals_adj[markers1$cluster == ct & markers1$pvals_adj > 0])))
  pval_limits <- 10^-c(0:abs(min_pval))
  
  po <- do.call(rbind, lapply(pval_limits, function(pl){
    sig_genes1 <- markers1$genes[markers1$cluster == ct & markers1$pvals_adj < pl ]
    sig_genes2 <- markers2$genes[markers2$cluster == ct & markers2$pvals_adj < pl ]
    
    intersect <- length(intersect(sig_genes1, sig_genes2))
    union <- length(union(sig_genes1, sig_genes2))
    
    if (intersect == 0 || union == 0){
      0
    }else{
      intersect / union
    }
    data.frame(cell_type = ct, pval_limit = pl, prop_overlap = intersect/ union, context = c('pancreas', 'whole mouse'), differential_genes = c(length(sig_genes1), length(sig_genes2)))
  }))
}))

# Plot the progressive change in intersect proportion for cell types over the
# list of p value thresholds. Plot the number differential genes on a second
# axis to provide context.

ggplot(prop_overlaps, aes(x = log10(pval_limit))) + 
  geom_line(aes(y = prop_overlap)) +
  geom_line(aes(y = log10(differential_genes)/5, color = context)) +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Proportion overlap significant marker genes (intersect/ union)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*5, name="Log10(Number of marker genes)")
  ) +
  theme_bw() + 
  xlab("log10(Adjusted p value limit)") + 
  ylab("Proportion overlap significant marker genes (intersect/ union)") +
  facet_wrap( ~ cell_type, scales = "free_x") + 
  labs(
    title = "Comparison of marker gene sets for pacreas cell types in Tabula Muris, pancreas context vs whole-mouse context",
    subtitle = "At progressively more relaxed adjusted p value thresholds, note the variable x axes.",
    caption = "Notes:\n\n1) Some cell types have few differential genes even at high p value thresholds, e.g. acinar cells\n2) Proportion of intersect tends towards 1 as the p value threshold is relaxed.\n3) Many more marker genes result from a whole mouse context than the pancreas-specific context\n4) Some cell types have non-trivial intersects across the p value range (e.g. type B pancreatic cell) implying there are context-independent marker genes, but since even at high stringency the\n    overall intersect is low (< 0.2), this set could not be identified from either marker set alone.\n5) Some cell types only develop significant intersect at low stringency (e.g. pancreatic PP cell), indicating that the marker genes are more context specific."
  ) + 
  theme(plot.caption = element_text(hjust = 0))
