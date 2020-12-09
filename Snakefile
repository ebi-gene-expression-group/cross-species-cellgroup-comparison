IN_DIR='inputs'
OUT_DIR='output'
OBO_FILE="%s/uberon.obo" % IN_DIR

EXPS=[config["exp1"], config["exp2"]]

rule all:
    input:
        #metas=expand("%s/meta/{exp}.tsv" % OUT_DIR, exp=EXPS)
        #cellgroup_mappings="%s/E-MTAB-5061_vs_E-ENAD-15.txt" % OUT_DIR
        dynamic(expand("%s/adata_subs/{exp}.{{organism_part}}.h5ad" % OUT_DIR, exp=EXPS))

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
        compare_terms.R {input.oboFile} {input.metas[0]} {input.metas[1]} organism_part_ontology organism_part {OUT_DIR}/meta_subs
        """

rule subset_data_toparts:
    conda:
         'envs/scanpy-scripts.yml'
    
    input:
        adata="%s/{exp}.project.h5ad" % IN_DIR,
        sub="{outdir}/meta_subs/{exp}.{organism_part}.tsv"        

    output:
        adata_sub="{outdir}/adata_subs/{exp}.{organism_part}.h5ad"

    shell:
        """
        scripts/subset_anndata.py {input.adata} {input.sub} {output.adata_sub}
        """
