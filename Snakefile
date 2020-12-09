IN_DIR='inputs'
OUT_DIR='output'
OBO_FILE="%s/uberon.obo" % IN_DIR

EXPS=[config["exp1"], config["exp2"]]
CELL_TYPE_FIELD=config['cell_type_field']

rule all:
    input:
        #metas=expand("%s/meta/{exp}.tsv" % OUT_DIR, exp=EXPS)
        #cellgroup_mappings="%s/E-MTAB-5061_vs_E-ENAD-15.txt" % OUT_DIR
        #dynamic(expand("%s/markers/{exp}.{{organism_part}}.tsv" % OUT_DIR, exp=EXPS))
        dynamic(expand("%s/markers/{exp}.{{organism_part}}.markers.h5ad" % OUT_DIR, exp=EXPS))

rule extract_metadata:
    conda:
         'envs/scanpy-scripts.yml'

    input:
        anndata = "%s/{exp}.project.h5ad" % IN_DIR

    output:
        meta = "{outdir}/meta/{exp}.tsv"

    shell:
        """
        scripts/extract_meta.py {input.anndata} {output.meta}        
        """

rule intersect_metadatas:
    conda:
         'envs/ontology_index.yml'
    
    input:
        oboFile = OBO_FILE,
        metas = expand("{{outdir}}/meta/{exp}.tsv", exp = EXPS)
    
    output:
        dynamic(expand("{{outdir}}/meta_subs/{exp}.{{organism_part}}.tsv", exp=EXPS))

    shell:
        """
        bin/compare_terms.R {input.oboFile} {input.metas[0]} {input.metas[1]} organism_part_ontology organism_part {OUT_DIR}/meta_subs
        """

rule subset_data_toparts:
    conda:
         'envs/scanpy-scripts.yml'
    
    input:
        adata="%s/{exp}.project.h5ad" % IN_DIR,
        sub="{outdir}/meta_subs/{exp}.{organism_part}.tsv"        

    output:
        adata_sub=temp("{outdir}/adata_subs/{exp}.{organism_part}.h5ad")

    shell:
        """
        scripts/subset_anndata.py {input.adata} {input.sub} {output.adata_sub}
        """

rule filter_cells:
    conda:
         'envs/scanpy-scripts.yml'
    
    input:
        adata_sub="{outdir}/adata_subs/{exp}.{organism_part}.h5ad"

    output:
        adata=temp("{outdir}/markers/{exp}.{organism_part}.cellfiltered.h5ad")

    shell:
        """
        scanpy-filter-cells --gene-name 'gene_symbols' --param 'c:n_genes' 400.0 \
            1000000000.0 --param 'c:n_counts' 0.0 1000000000.0  --input-format 'anndata' \
            {input.adata_sub}  --show-obj stdout --output-format anndata {output.adata}
        """
        
rule filter_genes:
    conda:
         'envs/scanpy-scripts.yml'
    
    input:
        adata="{outdir}/markers/{exp}.{organism_part}.cellfiltered.h5ad"

    output:
        adata=temp("{outdir}/markers/{exp}.{organism_part}.genefiltered.h5ad")

    shell:
        """
        scanpy-filter-genes --param 'g:n_cells' 3.0 1000000000.0 --input-format 'anndata' \
            {input.adata}  --show-obj stdout --output-format anndata {output.adata}
        """


rule normalise:
    conda:
         'envs/scanpy-scripts.yml'
    
    input:
        adata="{outdir}/markers/{exp}.{organism_part}.genefiltered.h5ad"

    output:
        adata=temp("{outdir}/markers/{exp}.{organism_part}.normalised.h5ad")

    shell:
        """
        scanpy-normalise-data --normalize-to '10000.0' --save-raw yes \
            --input-format 'anndata' {input.adata}  --show-obj stdout --output-format anndata \
            {output.adata}
        """

rule rgg:
    conda:
         'envs/scanpy-scripts.yml'
    
    input:
        adata="{outdir}/markers/{exp}.{organism_part}.normalised.h5ad"

    output:
        adata="{outdir}/markers/{exp}.{organism_part}.markers.h5ad",
        tsv="{outdir}/markers/{exp}.{organism_part}.markers.tsv"

    shell:
        """
        scanpy-find-markers --save "{output.tsv}" --n-genes '100' --groupby '{CELL_TYPE_FIELD}' \
            --key-added 'markers_{CELL_TYPE_FIELD}' --method 't-test_overestim_var' --use-raw \
             --reference 'rest' --filter-params 'min_in_group_fraction:0.0,max_out_group_fraction:1.0,min_fold_change:1.0' \
             --input-format 'anndata' {input.adata} --show-obj stdout \
            --output-format anndata {output.adata}
        """


