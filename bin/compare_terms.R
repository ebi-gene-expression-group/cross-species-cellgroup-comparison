#!/usr/bin/env Rscript

library(ontologyIndex)

cl <- commandArgs(trailingOnly = TRUE)
oboFile <- cl[1]
metaFile1 <- cl[2]
metaFile2 <- cl[3]
ontIdCol <- cl[4]
ontNameCol <- cl[5]
outPath <- cl[6]

# Take two sets of ontology annotations, determine if any terms in one list
# match, or are children of terms in another, and attempt to harmonise

ont <- ontologyIndex::get_OBO(oboFile, propagate_relationships = c('is_a', 'part_of'))

meta1 <- read.delim(metaFile1, stringsAsFactors = F)
meta2 <- read.delim(metaFile2, stringsAsFactors = F)

onts1 <- sub('_', ':', unique(basename(meta1[[ontIdCol]])))
onts2 <- sub('_', ':', unique(basename(meta2[[ontIdCol]])))

# 1. Get direct mappings

intersection <- intersect(onts1, onts2)
mappings = c()
mappings[intersection] = intersection

onts1 <- onts1[! onts1 %in% intersection]
onts2 <- onts2[! onts2 %in% intersection]

# 2. For the remaining terms, determine if they are found in the ancestor lists
# of the other dataset's terms

ancestors_onts_1 <- lapply(onts1, function(x) get_ancestors(ont, x))
ancestors_onts_2 <- lapply(onts2, function(x) get_ancestors(ont, x))

# For any term in the second list, map it to any ancestor term found in the
# first list

one_ancestors_of_two <- lapply(structure(onts1, names=onts1), function(x) onts2[unlist(lapply(ancestors_onts_2, function(y) x %in% y))])
one_ancestors_of_two <- one_ancestors_of_two[unlist(lapply(one_ancestors_of_two, function(x) length(x) > 0))]

for (ancestor in names(one_ancestors_of_two)){
  mappings[sub(':', '_', one_ancestors_of_two[[ancestor]])] <- sub(':', '_', ancestor)
  
  onts1 <- onts1[onts1 != ancestor]
  onts2 <- onts2[! onts2 %in% one_ancestors_of_two[[ancestor]]]
}

# For any remaining term in the first list, map it to any ancestor term found in
# the second list

two_ancestors_of_one <- lapply(structure(onts2, names = onts2), function(x) onts1[unlist(lapply(ancestors_onts_1, function(y) x %in% y))])
two_ancestors_of_one <- two_ancestors_of_one[unlist(lapply(two_ancestors_of_one, function(x) length(x) > 0))]

for (ancestor in names(two_ancestors_of_one)){
  mappings[sub(':', '_', two_ancestors_of_one[[ancestor]])] <- sub(':', '_', ancestor)
  
  onts2 <- onts2[onts2 != ancestor]
  onts1 <- onts1[! onts1 %in% two_ancestors_of_one[[ancestor]]]
}

# The mappings vector now describes the substitutions that must be made in the
# ontology terms across both lists

for (fromTerm in names(mappings)){
  toTerm <- mappings[fromTerm]
  toTermName <- ont$name[sub('_', ':', toTerm)]
  
  search_term <- paste0('/', names(mappings)[1])
  
  meta1_matches <- grep(search_term, meta1[[ontIdCol]])
  meta1[meta1_matches, ontIdCol] <- sub(search_term, paste0('/', mappings[fromTerm]), meta1[meta1_matches, ontIdCol])
  meta1[meta1_matches, ontNameCol] <- toTermName
  
  meta2_matches <- grep(search_term, meta2[[ontIdCol]])
  meta2[meta2_matches, ontIdCol] <- sub(search_term, paste0('/', mappings[fromTerm]), meta2[meta2_matches, ontIdCol])
  meta2[meta2_matches, ontNameCol] <- toTermName
}

# Now find the intersection

common_parts <- intersect(meta1$organism_part_ontology, meta2$organism_part_ontology)
meta1 <- meta1[meta1[[ontIdCol]] %in% common_parts,  ]
meta2 <- meta2[meta2[[ontIdCol]] %in% common_parts,  ]

# Output a metadata path 

metas1 <- split(meta1, meta1[[ontNameCol]])
metas2 <- split(meta2, meta2[[ontNameCol]])

dir.create(outPath, showWarnings=FALSE)

for (metaName in names(metas1)){
  write.table(metas1[[metaName]], file.path(outPath, paste0(basename(tools::file_path_sans_ext(metaFile1)), '.', gsub(' ', '_', metaName), '.tsv')), quote = FALSE, row.names = FALSE, sep = "\t")
  write.table(metas2[[metaName]], file.path(outPath, paste0(basename(tools::file_path_sans_ext(metaFile2)), '.', gsub(' ', '_', metaName), '.tsv')), quote = FALSE, row.names = FALSE, sep = "\t")
}




