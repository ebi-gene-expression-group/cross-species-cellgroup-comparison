nextflow.enable.dsl=2

// Extract metadata

process extract_metadata {

    publishDir 'outputs'

    conda 'envs/scanpy.yml'

    input:
      tuple val(id), path(anndata)
    output:
      tuple val(id), path("${id}.tsv")
    
    """
    #!/usr/bin/env python
    import scanpy as sc
 
    ad = sc.read("$anndata") 
    ad.obs.to_csv("${id}.tsv", sep="\t")
    """
}


// Relate and intersect organism parts, and split outby organism part

process intersect_metadatas {

    publishDir 'outputs'
    
    conda 'envs/ontology_index.yml'
    
    input:
        tuple path(metaFile1), path(metaFile2)
    output:
        tuple val("${metaFile1.baseName}"), path("intersected/1/*"), emit: intersections_1
        tuple val("${metaFile2.baseName}"), path("intersected/2/*"), emit: intersections_2

    """
    compare_terms.R ${params.oboFile} $metaFile1 $metaFile2 organism_part_ontology organism_part intersected
    """

}

process extract_orgpart {

    input:
        tuple val(expName), path(metaFile)
    output:
        tuple val(expName), stdout, path(metaFile)


    """
    echo -n "$metaFile" | awk -F'.' '{print \$2}' | tr -d '\n'
    """
}

process subset_to_parts {

    publishDir 'outputs'
    
    conda 'envs/scanpy.yml'
    
    input:
        tuple val(expName), val(orgPart), path(anndata), path(metaFile)
    output:
        tuple val(expName), val(orgPart), path("${expName}.${orgPart}.h5ad")

    """
    #!/usr/bin/env python
    import scanpy as sc
    import pandas as pd 

    ad = sc.read("$anndata") 
    obs = pd.read_csv('$metaFile', sep = "\\t", index_col = 0)
    ad[obs.index].write("${expName}.${orgPart}.h5ad")
    """        
}

// Re-do marker detection with scanpy-scripts

process rank_genes_groups {

    publishDir 'outputs/rergg', mode: 'copy'
    
    conda 'envs/scanpy-scripts.yml'
    
    input:
        tuple val(expName), val(orgPart), path(annData)

    output:
        tuple val(expName), val(orgPart), path("${expName}.${orgPart}.rgg.h5ad"), path("${expName}.${orgPart}.rgg.tsv")

    """

    scanpy-filter-cells --gene-name 'gene_symbols' --param 'c:n_genes' 400.0 \
        1000000000.0 --param 'c:n_counts' 0.0 1000000000.0  --input-format 'anndata' \
        $annData   --show-obj stdout --output-format anndata cellfiltered.h5ad

    scanpy-filter-genes --param 'g:n_cells' 3.0 1000000000.0 --input-format 'anndata' \
        cellfiltered.h5ad --show-obj stdout --output-format anndata genefiltered.h5ad

    # This probably isn't doing anything with use_raw below, but we'll do it just in case

    scanpy-normalise-data --normalize-to '10000.0' --save-raw yes \
        --input-format 'anndata' genefiltered.h5ad  --show-obj stdout --output-format anndata \
        normalised.h5ad

    scanpy-find-markers --save "${expName}.${orgPart}.rgg.tsv" --n-genes '100' --groupby '${params.cell_type_field}' \
        --key-added 'markers_${params.cell_type_field}' --method 't-test_overestim_var' --use-raw \
         --reference 'rest' --filter-params 'min_in_group_fraction:0.0,max_out_group_fraction:1.0,min_fold_change:1.0' \
         --input-format 'anndata' normalised.h5ad --show-obj stdout \
        --output-format anndata "${expName}.${orgPart}.rgg.h5ad" 
    """

}

process compare_experiments {

    input:
        tuple val(tissue), val(expId1), val(expId2), path(markers1), path(markers2)

    output:

    """
    compare_experiments.R ${markers1} ${markers2}
    """

}

// Compare marker gene lists to generate a list of potential relationships

workflow {

    exps = Channel.from(["${params.expid1}"], ["${params.expid2}"])
    ad1 = ["${params.expid1}", file("${params.anndata1}")]
    ad2 = ["${params.expid2}", file("${params.anndata2}")]

    meta1 = extract_metadata(ad1)
    meta2 = extract_metadata(ad2)


//Channel.from( ["${params.expid1}", file("${params.anndata1}")], ["${params.expid2}", file("${params.anndata2}")] )
    //species = Channel.from( ["${params.expid1}", "${params.species1}"], ["${params.expid2}", "${params.species2}"] )


    //metadatas = exps.combine(extract_metadata(ads), by: 0)
    //intersect_metadatas( metadatas.toList().map{r -> r.flatten()} ) 
    
    //exps.view()
    //metadatas.view()    

//metadatas.toList().map{r -> r.flatten()}.view()


//named_anndatas = extract_name(ads)
    //named_anndatas.view()

   // Channel.from(tuple("${params.species1}", file("${params.anndata1}"))).concat(Channel.from(tuple("${params.species2", file("${params.anndata2}")))).view()


    //ad1 = Channel.fromPath(params.anndata1).first()  
    //ad2 = Channel.fromPath(params.anndata2).first()
    
    //ad1 = Channel.from(["${params.species1}", file("${params.anndata1}"]))
    //ad2 = Channel.from(["${params.species2}", file("${params.anndata2}"]))
    //ads = ad1.concat(ad2)
    
    //ads.view()


    //species = Channel.from([["${params.species1}", ad1], [ "${params.species2}", ad2]]).transpose()

    //ads.view()
    //species.view()

    //ads.merge(species).view()


    //inputs = Channel.from([params.species1, ad1], [params.species2, ad2]])
    
    //ads = ad1.concat(ad2)

    //metadatas = extract_metadata(named_anndatas)

    //labelled_ad_subsets = extract_orgpart(intersect_metadatas.out.intersections_1.concat(intersect_metadatas.out.intersections_2)) 

    //named_anndatas.cross(labelled_ad_subsets).map{ r -> r.flatten()[0,3,1,4]}.view()

    //subset_to_parts(named_anndatas.cross(labelled_ad_subsets).map{ r -> r.flatten()[0,3,1,4]})
    //rgg = rank_genes_groups(subset_to_parts.out).groupTuple(by: 1).map{ r -> tuple(r[1], r[0][0], r[0][1], r[3][0], r[3][1] }


}
