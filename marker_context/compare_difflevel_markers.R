#!/usr/bin/env Rscript

library(ggplot2)

# Take two marker gene lists derived from different contexts, and see how the
# overlap in genes changes as we make p value more and more stringent.

cl <- commandArgs(trailingOnly = TRUE)

marker_file_1 <- cl[1]
context_1 <- cl[2]
marker_file_2 <- cl[3]
context_2 <- cl[4]
mode <- cl[5]
outfile <- cl[6]

#marker_file_1 <- "~/output/E-ENAD-15.pancreas.markers.tsv"
#marker_file_2 <- "~/output/E-ENAD-15.markers.tsv"
#context_1 <- 'pancreas'
#context_2 <- 'whole mouse'
#mode = 'pval'
#rank_coeff <- 3
#outfile <- '~/Desktop/marker_context_rank.png'

markers1 <- read.delim(marker_file_1, stringsAsFactors = FALSE)
markers2 <- read.delim(marker_file_2, stringsAsFactors = FALSE)

common_celltypes <- unique(intersect(markers1$cluster, markers2$cluster))

# Find the proporition of overlapping marker genes at progressivlely relaxed
# thrsholds, starting at the minimum adjusted p value for each cell type

prop_overlaps <- do.call(rbind, lapply(common_celltypes, function(ct){
  
  print(ct)
  
  if (mode == 'pval'){
    
    min_pval <- floor(log10(min(markers1$pvals_adj[markers1$cluster == ct & markers1$pvals_adj > 0])))
    thresholds <- 10^-c(0:abs(min_pval))
    
    do.call(rbind, lapply(thresholds, function(t){
      sig_genes1 <- markers1$genes[markers1$cluster == ct & markers1$pvals_adj < t ]
      sig_genes2 <- markers2$genes[markers2$cluster == ct & markers2$pvals_adj < t ]
      
      intersect <- length(intersect(sig_genes1, sig_genes2))
      union <- length(union(sig_genes1, sig_genes2))
      
      res <- data.frame(cell_type = ct, intersect_union = intersect/ union, intersect_denominator = intersect/ length(sig_genes1), context = c(context_1, context_2), differential_genes = c(length(sig_genes1), length(sig_genes2)))
      res[[mode]] <- t
      res
    }))
  }else{
    
    max_rank <- min(c(1000, min(max(markers1$rank[markers1$cluster == ct]), max(markers1$rank[markers1$cluster == ct]))))
    thresholds <- 0:max_rank
    
    po <- do.call(rbind, lapply(thresholds, function(t){
      sig_genes1 <- markers1$genes[markers1$cluster == ct & markers1$rank <  t ]
      do.call(rbind, lapply(c(1, 3, 5, 10), function(rank_coeff){
        sig_genes2 <- markers2$genes[markers2$cluster == ct & markers2$rank < min((t * rank_coeff), max_rank) ] 
        
        intersect <- length(intersect(sig_genes1, sig_genes2))
        union <- length(union(sig_genes1, sig_genes2))
        
        res <- data.frame(cell_type = ct, intersect_union = intersect/ union, intersect_denominator = intersect/ length(sig_genes1), context = c(context_1, context_2), differential_genes = c(length(sig_genes1), length(sig_genes2)))
        res$rank_coeff <- rank_coeff
        res[[mode]] <- t
        res
      }))
    }))
    po$rank_coeff <- factor(po$rank_coeff)
    po
  }
}))

colnames(prop_overlaps)[colnames(prop_overlaps) == 'intersect_denominator'] <- paste0('intersect_', context_1)

# Plot the progressive change in intersect proportion for cell types over the
# list of p value thresholds. Plot the number differential genes on a second
# axis to provide context.

if (mode == 'pval'){
  
  # Reshape to plot different intersects alongside
  
  plotdata <- reshape2::melt(prop_overlaps, id.var = c('cell_type', 'context', mode, 'differential_genes'))
  colnames(plotdata)[5:6] <- c('overlap_type', 'prop_overlap')
  
  p <- ggplot(plotdata, aes(x = log10(pval))) +
    geom_line(aes(y = prop_overlap, linetype = overlap_type)) +
    geom_line(aes(y = log10(differential_genes)/5, color = context)) +
    scale_y_continuous(
    
      # Features of the first axis
      name = "Proportion overlap significant marker genes",
    
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*5, name="Log10(Number of marker genes)")
    ) +
    xlab("log10(Adjusted p value limit)")
}else{
  
  # Reshape to plot different intersects alongside
  
  plotdata <- reshape2::melt(prop_overlaps, id.var = c('cell_type', 'context', mode, 'differential_genes', 'rank_coeff'))
  colnames(plotdata)[6:7] <- c('overlap_type', 'prop_overlap')
  
  p <- ggplot(plotdata, aes(x = log10(rank))) + 
    geom_line(aes(y = prop_overlap, linetype = overlap_type, color = rank_coeff)) +
    xlab("log10(marker gene rank)")
}

p <- p + 
  theme_bw() + 
  ylab("Proportion overlap significant marker gene") +
  facet_wrap( ~ cell_type, scales = "free_x") + 
  labs(
    title = paste0("Comparison of marker gene sets for ", context_1 ," cell types, ", context_1, " context vs ", context_2," context")
  ) + 
  theme(plot.caption = element_text(hjust = 0))

png(outfile, width = 1000, height = 800)
print(p)
dev.off()

