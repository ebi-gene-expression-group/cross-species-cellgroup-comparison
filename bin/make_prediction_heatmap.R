#!/usr/bin/env Rscript

cl <- commandArgs(trailingOnly = TRUE)

library(pheatmap)
library(grid)

prediction_file <- cl[1]
method_name <- cl[2]
metric_name <- cl[3]
outfile <- cl[4]

prediction <- read.delim(prediction_file, stringsAsFactors = FALSE, check.names = FALSE)[,1:3]
names <- colnames(prediction)[1:2]
colnames(prediction) <- c('id', 'variable', 'value')

prediction <- reshape2::acast(prediction, id ~ variable, value.var="value")
prediction <- prediction[sort(rownames(prediction)), sort(colnames(prediction))]
prediction[is.na(prediction)] <- 0

png(filename = outfile, width = 1000, height = 1000)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(
  prediction,
  main = paste0('Cell group match prediction, ', method_name , ' (', metric_name, ')'),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  angle_col = 315,
  fontsize = 15
)
setHook("grid.newpage", NULL, "replace")
grid.text(names[2], y=-0.07, gp=gpar(fontsize=16))
grid.text(names[1], x=-0.07, rot=90, gp=gpar(fontsize=16))
dev.off()




