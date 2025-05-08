rule contig_annotate__camper__annotate:
    """Annotate dereplicate genomes with CAMPER"""
    input:
        fa=CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.fa",
    output:
        annotation = CONTIG_CAMPER / "{assembly_id}" / "annotations.tsv",
        distillate = CONTIG_CAMPER / "{assembly_id}" / "distillate.tsv",
    log:
        CONTIG_CAMPER / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    container:
        docker["camper"]
    params:
        camper_dir = CONTIG_CAMPER,
        out_dir = lambda wildcards: f"{CONTIG_CAMPER}/{wildcards.assembly_id}",
        camper_db=features["databases"]["camper"],
        nvme=config["nvme_storage"],
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
        nvme = config["resources"]["nvme"]["small"]
    shell:
        """
        
        cp {params.camper_db}/* {params.nvme} \
        2>> {log} 1>&2
        
        tmpdir={params.out_dir}/tmp
        
        camper_annotate -i {input} \
                        -o $tmpdir \
                        --threads {threads} \
                        --camper_fa_db_loc {params.nvme}/CAMPER_blast.faa \
	                      --camper_fa_db_cutoffs_loc {params.nvme}/CAMPER_blast_scores.tsv \
	                      --camper_hmm_loc {params.nvme}/CAMPER.hmm  \
                        --camper_hmm_cutoffs_loc {params.nvme}/CAMPER_hmm_scores.tsv \
        2>> {log} 1>&2

        camper_distill  -a $tmpdir/annotations.tsv \
                        -o {output.distillate} \
	                      --camper_distillate {params.nvme}/CAMPER_distillate.tsv \
	      2>> {log} 1>&2
	      
	      # Move all tmp files one level up
        mv "$tmpdir"/* {params.out_dir}/

        """


rule contig_annotate__camper:
    """Run CAMPER on contig gene sequences."""
    input:
         [CONTIG_CAMPER / f"{assembly_id}" / "annotations.tsv" for assembly_id in ASSEMBLIES],
