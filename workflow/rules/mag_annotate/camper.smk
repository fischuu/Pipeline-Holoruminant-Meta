rule mag_annotate__camper__annotate:
    """Annotate dereplicate genomes with CAMPER"""
    input:
        dereplicated_genomes=DREP / "dereplicated_genomes.fa.gz",
        annotation=DRAM / "annotate" / "annotations.tsv",
    output:
        annotation=CAMPER / "annotations.tsv",
        distill=CAMPER / "distill.tsv"        
    log:
        annot = CAMPER / "annotate.log",
        distill = CAMPER / "distill.log",
    conda:
        "__environment__.yml"
    container:
        docker["camper"]
    params:
        camper_dir = CAMPER,
        camper_db=features["databases"]["camper"],
        nvme=config["nvme_storage"],
    threads: esc("cpus", "mag_annotate__camper__annotate")
    resources:
        runtime=esc("runtime", "mag_annotate__camper__annotate"),
        mem_mb=esc("mem_mb", "mag_annotate__camper__annotate"),
        cpus_per_task=esc("cpus", "mag_annotate__camper__annotate"),
        slurm_partition=esc("partition", "mag_annotate__camper__annotate"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'mag_annotate__camper__annotate')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("mag_annotate__camper__annotate"))
    shell:
        """
        cp {params.camper_db}/* {params.nvme} 2>> {log} 1>&2
        
        tmpdir={params.camper_dir}/tmp
        
        camper_annotate -i {input.dereplicated_genomes} \
                        -a {input.annotation} \
                        -o $tmpdir \
                        --threads {threads} \
                        --camper_fa_db_loc {params.nvme}/CAMPER_blast.faa \
	                      --camper_fa_db_cutoffs_loc {params.nvme}/CAMPER_blast_scores.tsv \
	                      --camper_hmm_loc {params.nvme}/CAMPER.hmm  \
                        --camper_hmm_cutoffs_loc {params.nvme}/CAMPER_hmm_scores.tsv \
                        2>> {log.annot} 1>&2
        
        camper_distill  -a $tmpdir/annotations.tsv \
                        -o {output.distill} \
	                      --camper_distillate {params.nvme}/CAMPER_distillate.tsv \
	      2>> {log} 1>&2
	      
	      # Move all tmp files one level up
        mv "$tmpdir"/* {params.camper_dir}/
        """

rule mag_annotate__camper:
    """Run CAMPER on dereplicated genomes."""
    input:
        rules.mag_annotate__camper__annotate.output,
