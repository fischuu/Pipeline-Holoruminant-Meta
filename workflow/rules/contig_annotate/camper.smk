rule contig_annotate__camper__annotate:
    """Annotate dereplicate genomes with CAMPER"""
    input:
        fa=CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.fa",
    output:
        annotation = CONTIG_CAMPER / "{assembly_id}" / "annotations.tsv",
        distillate = CONTIG_CAMPER / "{assembly_id}" / "distillate.tsv",
    log:
        CONTIG_CAMPER / "{assembly_id}.log",
    container:
        docker["camper"]
    params:
        camper_dir = CONTIG_CAMPER,
        out_dir = lambda wildcards: f"{CONTIG_CAMPER}/{wildcards.assembly_id}",
        camper_db=features["databases"]["camper"],
        nvme=config["nvme_storage"],
    threads: esc("cpus", "contig_annotate__camper__annotate")
    resources:
        runtime=esc("runtime", "contig_annotate__camper__annotate"),
        mem_mb=esc("mem_mb", "contig_annotate__camper__annotate"),
        cpus_per_task=esc("cpus", "contig_annotate__camper__annotate"),
        slurm_partition=esc("partition", "contig_annotate__camper__annotate"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "contig_annotate__camper__annotate", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("contig_annotate__camper__annotate"))
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
