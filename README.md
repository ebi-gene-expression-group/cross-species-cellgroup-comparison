# Comparing cell groupings between experiments across species

In Single-cell Expression Atlas we're interesting in relating cell groupings (clusters, cell types) between experiments and across species boundaries, which we can do via the 'marker' genes of each group. Since each set of marker genes is a product of the context in which it was derived (sorted cell population, sub-tissue, tissue, whole organism), that context must be matched for comparison of marker gene sets (via ortholog relationships) to be valid. With that in mind this workflow will:

 1. Take anndata objects from SCXA analysis for two experiments containing comparable 'organism' parts, even if they're not labelled at the same granularity. 
 2. Match the organism parts beween experiments, using the Uberon ontology to re-label where the granularity of organism part annotation is not consistent.
 3. Subset the experiments to only the common organism parts. 
 4. Re-filter, re-normalise, and re-derive marker genes using Scanpy.
 5. Derive mapped pairs of cell groupings between the two experiments.  
