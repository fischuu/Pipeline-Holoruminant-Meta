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
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'contig_annotate__camper__annotate')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("contig_annotate__camper__annotate"))
    shell:
        """
        
        # Previously we copied the camper database to a faster disc space, but as
        # it is rather small db, we maybe even can keep it in the original place
        # I'll test it to keep it there, in case the step is getting too I/O intense
        # or it takes too long to complete, we can copy the db again to a faster
        # disc space
        # cp {params.camper_db}/* {params.nvme} 2>> {log} 1>&2
        #
        # mkdir -p $tmpdir_camper
        # # Copy CAMPER DB locally to tmpdir_camper to avoid collisions
        # cp {params.camper_db}/CAMPER.hmm $tmpdir_camper/CAMPER.hmm
        # cp {params.camper_db}/CAMPER_hmm_scores.tsv $tmpdir_camper/CAMPER_hmm_scores.tsv
        # cp {params.camper_db}/CAMPER_blast.faa $tmpdir_camper/CAMPER_blast.faa
        # cp {params.camper_db}/CAMPER_blast_scores.tsv $tmpdir_camper/CAMPER_blast_scores.tsv
        # cp {params.camper_db}/CAMPER_distillate.tsv $tmpdir_camper/CAMPER_distillate.tsv

        tmpdir={params.out_dir}/tmp
        
        # THis is a strange camper requiredment, when it runs in snakemake and folders might be created...
        mkdir -p $tmpdir
        rm -rf $tmpdir
        
        camper_annotate -i {input} \
                        -o $tmpdir \
                        --threads {threads} \
                        --camper_fa_db_loc {params.camper_db}/CAMPER_blast.faa \
	                      --camper_fa_db_cutoffs_loc {params.camper_db}/CAMPER_blast_scores.tsv \
	                      --camper_hmm_loc {params.camper_db}/CAMPER.hmm  \
                        --camper_hmm_cutoffs_loc {params.camper_db}/CAMPER_hmm_scores.tsv \
        2>> {log} 1>&2

        camper_distill  -a $tmpdir/annotations.tsv \
                        -o {output.distillate} \
	                      --camper_distillate {params.camper_db}/CAMPER_distillate.tsv \
	      2>> {log} 1>&2
	      
	      # Move all tmp files one level up
        mv "$tmpdir"/* {params.out_dir}/

        """


rule contig_annotate__camper:
    """Run CAMPER on contig gene sequences."""
    input:
         [CONTIG_CAMPER / f"{assembly_id}" / "annotations.tsv" for assembly_id in ASSEMBLIES],
