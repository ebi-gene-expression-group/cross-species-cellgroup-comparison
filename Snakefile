import itertools
configfile: "config.yaml"

IN_DIR='inputs'
OUT_DIR='output'

wildcard_constraints:
    exp = "[^\.]+",
    exp1 = "[^\.]+",
    exp2 = "[^\.]+",
    organism_part = "[^\.]+",
    method = "[^\.]+",
    pval_limit = "[^\-]+",
    n_genes = "[^\-]+",
    min_overlap = "[^\-]+",
    samap_threshold = "[^\-]+",
    subdir = "[^\/]+"

pval_limits = [ str(x) for x in config.get('compare_experiments').get('pval_limit') ]
n_gene_limits = [ str(x) for x in config.get('compare_experiments').get('n_genes') ]
min_overlap_limits = [ str(x) for x in config.get('compare_experiments').get('min_overlap') ]
samap_thresholds = [ str(x) for x in config.get('compare_experiments').get('samap_threshold') ]

EXPS = [ config.get(exp).get('id') for exp in ('exp1', 'exp2') ]
MARKERS_PARAMSETS = [ '-'.join(x) for x in list(itertools.product(pval_limits, n_gene_limits, min_overlap_limits)) ]

IN_DIR='inputs'
OUT_DIR='output'
COMP_OUT_DIR='%s/%s_vs_%s' % (OUT_DIR, EXPS[0], EXPS[1])


rule all:
    input:
        "%s/.done" % COMP_OUT_DIR        


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

    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000

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
        celltypes = "{prefix}.celltypes.tsv"

    params:
        cell_type_field = config.get('cell_type_field')        
    
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    
    shell:
        """
        index=$(head -n 1 {input.meta} | tr '\\t' '\\n' | grep -n '^{params.cell_type_field}$' | awk -F':' '{{print $1}}')
        tail -n +2 {input.meta} | cut -f $index | sort | uniq > {output.celltypes}
        """

checkpoint intersect_metadatas:
    conda:
         'envs/ontology_index.yml'
    
    input:
        ontology_file="%s/%s" % (IN_DIR, config.get('ontology_file')),
        meta1 = '{out_dir}/{exp1}.meta.tsv',
        meta2 = '{out_dir}/{exp2}.meta.tsv'

    output:
        submetas = directory('{out_dir}/{exp1}_vs_{exp2}/submetas/')
    
    shell:
        """
        bin/compare_terms.R {input.ontology_file} {input.meta1} {input.meta2} \
            organism_part_ontology organism_part {output.submetas}
        """

rule subset_data_toparts:
    conda:
         'envs/scanpy-scripts.yml'

    priority: 1    

    input:
        adata="{outdir}/{prefix}.h5ad",
        sub="%s/submetas/{prefix}.{organism_part}.submeta.tsv" % COMP_OUT_DIR       

    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    
    output:
        adata_sub=temp("{outdir}/{prefix}.{organism_part}.sub.h5ad")

    shell:
        """
        scripts/subset_anndata.py "{input.adata}" "{input.sub}" "{output.adata_sub}"
        """

# When we want to make markers from the un-split h5ad

rule whole_sub:
    input:
        adata="{prefix}.sub.h5ad",

    output:
        adata_sub=temp("{prefix}.sub.h5ad")

    shell:
        """
        ln -s {input.adata} {input.adata_sub}
        """


rule filter_cells:
    conda:
         'envs/scanpy-scripts.yml'
    
    input:
        adata_sub="{prefix}.sub.h5ad"

    output:
        adata=temp("{prefix}.cellfiltered.h5ad")

    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    
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

    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    
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

    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    
    shell:
        """
        scanpy-normalise-data --normalize-to '10000.0' --save-raw yes \
            --input-format 'anndata' {input.adata}  --show-obj stdout --output-format anndata \
            {output.adata}
        """

rule find_markers:
    conda:
         'envs/scanpy-scripts.yml'
    
    input:
        adata="{prefix}.normalised.h5ad"

    output:
        adata="{prefix}.markers.h5ad",
        tsv="{prefix}.markers.tsv"

    params:
        cell_type_field = config.get('cell_type_field'),   
        n_genes = '' if config.get('n_genes') is None else "--n-genes '%s'" % config.get('n_genes') 

    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64000
    
    shell:
        """
        scanpy-find-markers --save "{output.tsv}" --groupby \
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
        comp="{prefix}/{exp1}_vs_{exp2}.{organism_part}/markeroverlap.{pval_limit}-{n_genes}.prediction.tsv"     
    
    params:
        species1=config.get('exp1').get('species'),
        species2=config.get('exp2').get('species'),

    shell:
        """
        bin/compare_experiments.R {input.exp1} {params.species1} {input.exp2} \
            {params.species2} {input.ortholog_mapping_file} {wildcards.pval_limit} \
            {wildcards.n_genes} {output.comp}
        """

rule make_prediction_heatmap:
    conda:   
      'envs/pheatmap.yml'
    
    input:
        comp="{outdir}/{subdir}/{method}.{predict_params}.prediction.tsv",
        
    output:
        heatmap="{outdir}/{subdir}/{subdir}.{method}.{predict_params}.prediction.heatmap.png",
    
    shell:
        """
        bin/make_prediction_heatmap.R {input.comp} {wildcards.method} '{wildcards.method} score' {output.heatmap}
        """    


rule filter_prediction:
    input:
        comp="{prefix}/{method}.{params}.prediction.tsv"

    output:
        filtered="{prefix}/{method}.{params}.prediction.filtered_{threshold}.tsv"

    shell:
        """
        awk -F"\\t" -v threshold={wildcards.threshold} '{{ if($3 >= threshold) {{ print }}}}' {input.comp} > {output.filtered}
        if [ ! -s {output.filtered} ]; then
            echo "Did not find any cell group matches passing threshold of {wildcards.threshold} for {wildcards.method}" 1>&2
            exit 1
        fi
        """

rule compare_celltypes:
    input:
        exp1 = "{outdir}/{exp1}.{organism_part}.sub.celltypes.tsv",
        exp2 = "{outdir}/{exp2}.{organism_part}.sub.celltypes.tsv"

    output:
        txt=temp("{outdir}/{exp1}_vs_{exp2}.{organism_part}/celltypecomp.txt")     
    
    shell:
        """
        echo -e "## {wildcards.exp1} {wildcards.organism_part} cell types:\n" > {output.txt}
        cat {input.exp1} | sed 's/$/  /' | sed 's/^/ - /g' >> {output.txt}
        echo -e "\n" >> {output.txt}

        echo -e "## {wildcards.exp2} {wildcards.organism_part} cell types:\n" >> {output.txt}
        cat {input.exp2} | sed 's/$/  /' | sed 's/^/ - /g' >> {output.txt}
        echo -e "\n" >> {output.txt}

        echo -e "## Common cell types:\n" >> {output.txt}
        comm -12 {input.exp1} {input.exp2} | sed 's/$/  /' | sed 's/^/ - /g' >> {output.txt}
        """

rule compare_celltypes_prediction:
    input:
        exp1 = "{outdir}/{exp1}.{organism_part}.sub.celltypes.tsv",
        exp2 = "{outdir}/{exp2}.{organism_part}.sub.celltypes.tsv",
        heatmap="{outdir}/{exp1}_vs_{exp2}.{organism_part}/{exp1}_vs_{exp2}.{organism_part}.{method}.{predict_params}.prediction.heatmap.png",
        comp="{outdir}/{exp1}_vs_{exp2}.{organism_part}/{method}.{predict_params}.prediction.filtered_{threshold}.tsv",

    output:
        md=temp("{outdir}/{exp1}_vs_{exp2}.{organism_part}/{method}.{predict_params}.filtered_{threshold}.predictcomp.md")     
    
    shell:
        """
        comm -12 {input.exp1} {input.exp2} > intersect.txt
        intersect=$(cat intersect.txt | wc -l)
        correct=$(comm -12 intersect.txt <(cat {input.comp} | awk -F'\\t' '!_[$1]++' | awk -F '\t' '{{if($1==$2) print $0}}' | awk -F'\\t' '{{print $1}}' | sort | uniq) | sed 's/$/  /')
        correct_count=$(echo -e "$correct" | wc -l)
        echo -e "![{wildcards.method} prediction heatmap]($(basename {input.heatmap}))" > {output.md}
        echo -e "$correct_count of $intersect known intersecting cell types were predicted as top match by marker gene composition:  \n\n$correct\n" >> {output.md}
        almost_correct=$(comm -12 intersect.txt <(cat {input.comp} | awk -F '\t' '{{if($1==$2) print $0}}' | awk -F'\\t' '{{print $1}}' | sort | uniq))
        almost_correct_count=$(echo -e "$almost_correct" | wc -l)
        echo -e "$almost_correct_count of $intersect known intersecting cell types were predicted as a match (at any rank).  \n" >> {output.md}
        """

rule samap_config:
    conda:
         'envs/yaml.yml'
    
    input:
        adata1="{dir}/{exp1}.{organism_part}.sub.h5ad",
        adata2="{dir}/{exp2}.{organism_part}.sub.h5ad"
   
    output:
        yaml=temp("{dir}/{exp1}_vs_{exp2}.{organism_part}.samap_config.yaml")

    script:
        "scripts/make_samap_config.py" 

rule run_samap:
    input:
        adata1="{outdir}/{exp1}.{organism_part}.sub.h5ad",
        adata2="{outdir}/{exp2}.{organism_part}.sub.h5ad",
        transcriptome1=config.get('exp1').get('transcriptome'),
        transcriptome2=config.get('exp2').get('transcriptome'),
        config="{outdir}/{exp1}_vs_{exp2}.{organism_part}.samap_config.yaml"

    output:
        cell_type_map="{outdir}/{exp1}_vs_{exp2}.{organism_part}/samap_celltype_map.tsv",
        cell_type_map_png="{outdir}/{exp1}_vs_{exp2}.{organism_part}/samap_celltype_map.png"

    shell:
        """
        snakemake -s {workflow.basedir}/samap-workflow/workflow/Snakefile --profile lsf --rerun-incomplete --config configfile={input.config}
        mv {wildcards.outdir}/{wildcards.exp1}_vs_{wildcards.exp2}.{wildcards.organism_part}/*celltype_map.tsv {output.cell_type_map}
        mv {wildcards.outdir}/{wildcards.exp1}_vs_{wildcards.exp2}.{wildcards.organism_part}/*celltype_map_heatmap.png {output.cell_type_map_png} 
        """

# Get SAMap predictions in format we can use 

rule parse_samap_predictions:
    conda:
         'envs/reshape2.yml'
    
    input:
        cell_type_map="{outdir}/{exp1}_vs_{exp2}.{organism_part}/samap_celltype_map.tsv",

    output:
        comp="{outdir}/{exp1}_vs_{exp2}.{organism_part}/samap.defaults.prediction.tsv"     

    params:
        species1=config.get('exp1').get('species'),
        species2=config.get('exp2').get('species'),
    
    shell:
        """
        {workflow.basedir}/bin/reformat_scmap_table.R {input.cell_type_map} {params.species1} {params.species2} {output.comp}
        """
        

rule report_comparison:
    input:
        comp_markers="{outdir}/{exp1}_vs_{exp2}.{organism_part}/markeroverlap.{pval_limit}-{n_genes}.prediction.filtered_{min_overlap}.tsv",
        predictcomp_markers="{outdir}/{exp1}_vs_{exp2}.{organism_part}/markeroverlap.{pval_limit}-{n_genes}.filtered_{min_overlap}.predictcomp.md",
        comp_samap="{outdir}/{exp1}_vs_{exp2}.{organism_part}/samap.defaults.prediction.filtered_{samap_threshold}.tsv",
        predictcomp_samap="{outdir}/{exp1}_vs_{exp2}.{organism_part}/samap.defaults.filtered_{samap_threshold}.predictcomp.md",
        celltypes="{outdir}/{exp1}_vs_{exp2}.{organism_part}/celltypecomp.txt",
        samap_result="{outdir}/{exp1}_vs_{exp2}.{organism_part}/samap_celltype_map.tsv" 

    output:
        report="{outdir}/{exp1}_vs_{exp2}.{organism_part}/markers-{pval_limit}-{n_genes}-{min_overlap}.samap-{samap_threshold}.report.md"
    shell:
        """
        echo -e "# Known composition of inputs\n\n" > {output.report}
        cat {input.celltypes} >> {output.report}

        echo -e "\n# Cell group matches based on marker genes:\n" >> {output.report}
        echo -e "\n## Parameters  \n\n - Maximum p value: {wildcards.pval_limit}  \n - Minimum proportion overlap: {wildcards.min_overlap}  \n" >> {output.report}
        echo -e "## Results \n" >> {output.report}
        cat {input.predictcomp_markers} >> {output.report}
        cat {input.comp_markers} | sed 's/\t/ | /g' | sed 's/^/| /g' | sed 's/$/ |  /g' > tab.tmp
        head -n 1 tab.tmp >> {output.report}
        echo -e "| --- | --- | --- | --- | --- | --- |" >> {output.report}
        tail -n +2 tab.tmp >> {output.report}
        rm tab.tmp

        echo -e "\n# Cell group matches based on SAMap results:\n" >> {output.report}
        echo -e "\n## Parameters  \n\n - SAMap minimum score threshold: {wildcards.samap_threshold}  \n" >> {output.report}
        echo -e "## Results \n" >> {output.report}
        cat {input.predictcomp_samap} >> {output.report}
        cat {input.comp_samap} | sed 's/\t/ | /g' | sed 's/^/| /g' | sed 's/$/ |  /g' > tab.tmp
        head -n 1 tab.tmp >> {output.report}
        echo -e "| --- | --- | --- |" >> {output.report}
        tail -n +2 tab.tmp >> {output.report}
        rm tab.tmp


        """

# Use the checkpointed intersect_metadatas step to work out what common
# organism parts we had and generate a list of outputs

def organism_part_outputs(wildcards):
    checkpoint_output = checkpoints.intersect_metadatas.get(**wildcards).output[0]
    organism_parts, = glob_wildcards(os.path.join(checkpoint_output, "%s.{organism_part}.submeta.tsv" % EXPS[0]))
    return expand("%s/%s_vs_%s.{organism_part}/markers-{pval_limit}-{n_genes}-{min_overlap}.samap-{samap_threshold}.report.md" % (OUT_DIR, config.get('exp1').get('id'), config.get('exp2').get('id')),
           pval_limit=pval_limits,
           n_genes = n_gene_limits,
           min_overlap = min_overlap_limits,
           samap_threshold = samap_thresholds,
           organism_part = organism_parts)

rule aggregate_reports:
    input:
        organism_part_outputs
    output:
        done="{out_dir}/{exp1}_vs_{exp2}/.done"

    shell:
        """
        touch {output.done}
        """    
