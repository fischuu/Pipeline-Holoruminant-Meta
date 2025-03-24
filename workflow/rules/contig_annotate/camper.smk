rule contig_annotate__camper__annotate:
    """Annotate dereplicate genomes with CAMPER"""
    input:
        fa=CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.fa",
    output:
        annotation = CONTIG_CAMPER / "{assembly_id}"
    log:
        CONTIG_CAMPER / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    container:
        docker["camper"]
    params:
        out_dir=CAMPER,
        camper_db=features["databases"]["camper"],
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        camper_annotate -i {input} \
                        -o {output} \
                        --threads {threads} \
                        --camper_fa_db_loc {params.camper_db}/CAMPER_blast.faa \
	                      --camper_fa_db_cutoffs_loc {params.camper_db}/CAMPER_blast_scores.tsv \
	                      --camper_hmm_loc {params.camper_db}/CAMPER.hmm  \
                        --camper_hmm_cutoffs_loc {params.camper_db}/CAMPER_hmm_scores.tsv \
                        camper_annotate -i my_genes.faa -o my_output \
	
        2>> {log} 1>&2

        #camper_distill  -a <path to annotations.tsv> -o <name of output.tsv>
        #camper_distill  -a my_output/annotations.tsv -o my_output/distillate.tsv \
	      #--camper_distillate CAMPER_distillate.tsv
        """

rule contig_annotate__camper:
    """Run CAMPER on contig gene sequences."""
    input:
         [CONTIG_CAMPER / f"{assembly_id}" for assembly_id in ASSEMBLIES],
