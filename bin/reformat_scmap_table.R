#!/usr/bin/env Rscript

# Reformat the square samap table to a long one we can process similarly to the
# marker-derived one

cl <- commandArgs(trailingOnly = TRUE)

samap_file <- cl[1]
species1 <- cl[2]
species2 <- cl[3]
outfile <- cl[4]

samap <- reshape2::melt(read.delim(samap_file))
samap$variable <- sub('[^_]+_(.*)', '\\1', samap$variable)
samap$X <- sub('[^_]+_(.*)', '\\1', samap$X)

samap <- samap[order(samap$value, decreasing = TRUE), ]
samap <- samap[,c(2,1,3)]
samap$variable <- gsub('\\.', ' ', samap$variable)
colnames(samap) <- c(paste(c(species1, species2), 'cluster'), 'score')
write.table(samap, file = outfile, quote = FALSE, sep="\t", row.names = FALSE)
