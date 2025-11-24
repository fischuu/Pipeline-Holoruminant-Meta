rule mag_annotate__camper__annotate:
    """Annotate dereplicate genomes with CAMPER"""
    input:
        genes=DRAM / "annotate" / "genes.faa",
        annotation=DRAM / "annotate" / "annotations.tsv",
    output:
        annotation=CAMPER / "annotations.tsv",
        distill=CAMPER / "distill.tsv"        
    log:
        CAMPER / "camper.log",
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
        # Previously we copied the camper database to a faster disc space, but as
        # it is rather small db, we maybe even can keep it in the original place
        # I'll test it to keep it there, in case the step is getting too I/O intense
        # or it takes too long to complete, we can copy the db again to a faster
        # disc space
        # cp {params.camper_db}/* {params.nvme} 2>> {log} 1>&2
        
        tmpdir={params.camper_dir}/tmp
        
        # THis is a strange camper requiredment, when it runs in snakemake and folders might be created...
        mkdir -p $tmpdir
        rm -rf $tmpdir

        camper_annotate -i {input.genes} \
                        -o $tmpdir \
                        --threads {threads} \
                        --camper_fa_db_loc {params.camper_db}/CAMPER_blast.faa \
	                      --camper_fa_db_cutoffs_loc {params.camper_db}/CAMPER_blast_scores.tsv \
	                      --camper_hmm_loc {params.camper_db}/CAMPER.hmm  \
                        --camper_hmm_cutoffs_loc {params.camper_db}/CAMPER_hmm_scores.tsv \
                        2>> {log} 1>&2
        
        camper_distill  -a $tmpdir/annotations.tsv \
                        -o {output.distill} \
	                      --camper_distillate {params.camper_db}/CAMPER_distillate.tsv \
	      2>> {log} 1>&2
	      
	      # Move all tmp files one level up
        mv "$tmpdir"/* {params.camper_dir}/
        """

rule mag_annotate__camper:
    """Run CAMPER on dereplicated genomes."""
    input:
        rules.mag_annotate__camper__annotate.output,
