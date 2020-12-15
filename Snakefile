configfile: "config.yaml"

IN_DIR='inputs'
OUT_DIR='output'

EXPS=[config.get("exp1").get('id'), config.get("exp2").get("id")]

wildcard_constraints:
    exp = "[^\.]",

rule all:
    input:
        dynamic("%s/%s_vs_%s.{organism_part}.celltypecomp.txt" % (OUT_DIR, config.get('exp1').get('id'), config.get('exp2').get('id')))     
        #dynamic("%s/%s_vs_%s.{organism_part}.tsv" % (OUT_DIR, config.get('exp1').get('id'), config.get('exp2').get('id')))     
        #dynamic(expand("%s/{exp}.{{organism_part}}.markers.h5ad" % OUT_DIR, exp=EXPS))
        #dynamic(expand("%s/{exp}.{{organism_part}}.submeta.tsv" % OUT_DIR, exp=EXPS))

rule symlink_inputs:
    input:
        '%s/{file}.project.h5ad' % IN_DIR
    output:
        '%s/{file}.h5ad' % OUT_DIR
    shell:
        """
        cp $(pwd)/{input} $(pwd)/{output}
        """

rule extract_metadata:
    conda:
         'envs/scanpy-scripts.yml'

    input:
        anndata = "{prefix}.h5ad"

    output:
        meta = temp("{prefix}.meta.tsv")

    shell:
        """
        scripts/extract_meta.py {input.anndata} {output.meta}        
        """

rule extract_celltypes:
    conda:
         'envs/scanpy-scripts.yml'

    input:
        meta = "{prefix}.meta.tsv",

    output:
        celltypes = temp("{prefix}.celltypes.tsv")

    params:
        cell_type_field = config.get('cell_type_field')        
    
    shell:
        """
        index=$(head -n 1 {input.meta} | tr '\\t' '\\n' | grep -n '^{params.cell_type_field}$' | awk -F':' '{{print $1}}')
        tail -n +2 {input.meta} | cut -f $index | sort | uniq > {output.celltypes}
        """

rule intersect_metadatas:
    conda:
         'envs/ontology_index.yml'
    
    input:
        metas = expand("{{outdir}}/{exp}.meta.tsv", exp = EXPS),
        ontology_file="%s/%s" % (IN_DIR, config.get('ontology_file'))
    
    output:
        dynamic(expand("{{outdir}}/{exp}.{{organism_part}}.submeta.tsv", exp=EXPS))

    shell:
        """
        bin/compare_terms.R {input.ontology_file} {input.metas[0]} {input.metas[1]} \
            organism_part_ontology organism_part {OUT_DIR}
        """

rule subset_data_toparts:
    conda:
         'envs/scanpy-scripts.yml'

    priority: 1    

    input:
        adata="{prefix}.h5ad",
        sub="{prefix}.{organism_part}.submeta.tsv"        

    output:
        adata_sub=temp("{prefix}.{organism_part}.sub.h5ad")

    shell:
        """
        scripts/subset_anndata.py "{input.adata}" "{input.sub}" "{output.adata_sub}"
        """

rule filter_cells:
    conda:
         'envs/scanpy-scripts.yml'
    
    input:
        adata_sub="{prefix}.sub.h5ad"

    output:
        adata=temp("{prefix}.cellfiltered.h5ad")

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
        adata="{prefix}.cellfiltered.h5ad"

    output:
        adata=temp("{prefix}.genefiltered.h5ad")

    shell:
        """
        scanpy-filter-genes --param 'g:n_cells' 3.0 1000000000.0 --input-format 'anndata' \
            {input.adata}  --show-obj stdout --output-format anndata {output.adata}
        """


rule normalise:
    conda:
         'envs/scanpy-scripts.yml'
    
    input:
        adata="{prefix}.genefiltered.h5ad"

    output:
        adata=temp("{prefix}.normalised.h5ad")

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
        adata="{prefix}.normalised.h5ad"

    output:
        adata=temp("{prefix}.markers.h5ad"),
        tsv=temp("{prefix}.markers.tsv")

    params:
        cell_type_field = config.get('cell_type_field')        

    shell:
        """
        scanpy-find-markers --save "{output.tsv}" --n-genes '100' --groupby \
            '{params.cell_type_field}' --key-added 'markers_{params.cell_type_field}' --method \
            't-test_overestim_var' --use-raw --reference 'rest' --filter-params \
            'min_in_group_fraction:0.0,max_out_group_fraction:1.0,min_fold_change:1.0' \
            --input-format 'anndata' {input.adata} --show-obj stdout --output-format \
            anndata {output.adata}
        """

rule compare_experiments:
    input:
        exp1="{prefix}/{exp1}.{organism_part}.markers.tsv",
        exp2="{prefix}/{exp2}.{organism_part}.markers.tsv",
        ortholog_mapping_file="%s/%s" % (IN_DIR, config.get('ortholog_mapping_file'))

    output:
        comp="{prefix}/{exp1}_vs_{exp2}.{organism_part}.tsv"     
    
    params:
        species1=config.get('exp1').get('species'),
        species2=config.get('exp2').get('species'),
        min_overlap=config.get('compare_experiments').get('min_overlap'),
        pval_limit=config.get('compare_experiments').get('pval_limit')

    shell:
        """
        bin/compare_experiments.R {input.exp1} {params.species1} {input.exp2} \
            {params.species2} {input.ortholog_mapping_file} {params.pval_limit} \
            {params.min_overlap} {output.comp}
        """

rule compare_celltypes:
    input:
        exp1 = "{outdir}/{exp1}.{organism_part}.sub.celltypes.tsv",
        exp2 = "{outdir}/{exp2}.{organism_part}.sub.celltypes.tsv"

    output:
        txt="{outdir}/{exp1}_vs_{exp2}.{organism_part}.celltypecomp.txt"     
    
    shell:
        """
        echo -e "## {wildcards.exp1} {wildcards.organism_part} cell types:\n" > {output.txt}
        cat {input.exp1} >> {output.txt}
        echo -e "\n\n" >> {output.txt}

        echo -e "## {wildcards.exp1} {wildcards.organism_part} cell types:\n" >> {output.txt}
        cat {input.exp2} >> {output.txt}
        echo -e "\n\n" >> {output.txt}

        echo -e "Common cell types:\n" >> {output.txt}
        comm -12 {input.exp1} {input.exp2} >> {output.txt}
        """

rule report_comparison:
    input:
        comp="{outdir}/{exp1}_vs_{exp2}.{organism_part}.tsv",
        celltypes="{outdir}/{exp1}_vs_{exp2}.{organism_part}.celltypecomp.txt"

    output:
        report="{outdir}/{exp1}_vs_{exp2}.{organism_part}.report.md"

    shell:
        """
        cat {input.comp} > {output.report}
        echo -e "# Matches predicted from marker genes:" >> {output.report}
        cat {input.txt} >> {output.report}
        """






