# Best parameters for marker selection

What's the best way of selecting markers for cross-species comparison of cell groups? Should we just take the top 'n', all genes matching a p value threshold? I iterated over three parameters with some test values:

 - minimum cluster overlap to be considered a match: [ 0.05, 0.1, 0.2 ]
 - maximum adjusted p value: 1[ (no limit), 0.1, 0.05, 0.01, 0.001 ]
 - number of genes selected: [ 10, 100, 1000 ] 

For a total of 45 parameter combinations. I tried these parameters with the following comparisons:

 - E-GEOD-125970 vs E-ENAD-15 (colon)
 - E-GEOD-83139 vs E-ENAD-15 (pancreas)
 - E-HCAD-10 vs E-ENAD-15 (kidney)
 - E-HCAD-1 vs E-ENAD-15 (lung)
 - E-HCAD-1 vs E-ENAD-15 (spleen)
 - E-MTAB-5061 vs E-ENAD-15 (pancreas)
 - E-MTAB-8410 vs E-ENAD-15 (ascending colon)

... and determined which parameter sets produced the maximal correctly predicted interseting gene sets.


# Results

The parameter sets (max pval-n genes-overlap proportion) top-ranked in more than one comparison are:

 - 0.01-100-0.05 (2 comparisons)
 - 0.05-100-0.05 (2 comparisons)
 - 0.1-100-0.05 (2 comparisons)
 - 1-1000-0.2 (2 comparisons)
 - 1-1000-0.1 (4 comparisons)
 - 1-1000-0.05 (5 comparisons)
 - 1-100-0.05 (5 comparisons)

Clearly this is a small set of examples, but it seems like selecting the top 100 genes without a p value filter, and requring a 5 percent overlap is a sensible starting position. Coincidentally, this is currently how we select marker genes in our single-cell analyses.
