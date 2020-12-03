nextflow.enable.dsl=2

// Input both objects

// Extract first parat of name

process extract_name {

    input:
        path anndata
    output:
        tuple stdout, path(anndata)

    """
    echo -e $anndata | awk -F'.' '{print \$1}' | tr -d '\n'  
    """
}


// Extract metadata

process extract_metadata {

    publishDir 'outputs'

    conda 'envs/scanpy.yml'

    input:
      tuple val(id), path(anndata)
    output:
      path "${id}.tsv"
    
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
        tuple val("${metaFile1.simpleName}"), path("intersected/1/*"), emit: intersections_1
        tuple val("${metaFile2.simpleName}"), path("intersected/2/*"), emit: intersections_2

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
    scanpy-find-markers --save "${expName}.${orgPart}.rgg.tsv" --n-genes '100' --groupby '${params.cell_type_field}' \
        --key-added 'markers_${params.cell_type_field}' --method 't-test_overestim_var' --use-raw \
         --reference 'rest' --filter-params 'min_in_group_fraction:0.0,max_out_group_fraction:1.0,min_fold_change:1.0' \
         --input-format 'anndata' $annData --show-obj stdout --output-format anndata "${expName}.${orgPart}.rgg.h5ad" 
    """

}

// Compare marker gene lists to generate a list of potential relationships

workflow {
    ad1 = Channel.fromPath(params.anndata1).first()  
    ad2 = Channel.fromPath(params.anndata2).first()

    ads = ad1.concat(ad2)

    inputs = Channel.fromList([ad1, ad2])
    named_anndatas = extract_name(ads)
    metadatas = extract_metadata(named_anndatas)

    intersect_metadatas(metadatas.toList())    
    labelled_ad_subsets = extract_orgpart(intersect_metadatas.out.intersections_1.concat(intersect_metadatas.out.intersections_2)) 

    //named_anndatas.cross(labelled_ad_subsets).map{ r -> r.flatten()[0,3,1,4]}.view()

    subset_to_parts(named_anndatas.cross(labelled_ad_subsets).map{ r -> r.flatten()[0,3,1,4]})
    rank_genes_groups(subset_to_parts.out).groupTuple(by: 1).view()
}