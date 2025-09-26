rule mag_annotate__eggnog:
    """Run eggnog over the dereplicated mags"""
    input:
        contigs=DREP / "dereplicated_genomes.fa.gz",
    output:
        directory(EGGNOG) ,
    log:
        protected(EGGNOG / "eggnog.log"),
    container:
        docker["eggnog"]
    params:
        out_dir=EGGNOG,
        db=features["databases"]["eggnog"],
        prefix="eggnog"
    threads: esc("cpus", "mag_annotate__eggnog")
    resources:
        runtime=esc("runtime", "mag_annotate__eggnog"),
        mem_mb=esc("mem_mb", "mag_annotate__eggnog"),
        cpus_per_task=esc("cpus", "mag_annotate__eggnog"),
        slurm_partition=esc("partition", "mag_annotate__eggnog"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'mag_annotate__eggnog')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("mag_annotate__eggnog"))
    shell:
        """
         cp -r {params.db}/* $TMPDIR  2>> {log} 1>&2;
        
        emapper.py -m diamond \
                   --data_dir $TMPDIR \
                   --itype metagenome \
                   --genepred prodigal \
                   --dbmem \
                   --no_annot \
                   --no_file_comments \
                   --cpu {threads} \
                   -i {input.contigs} \
                   --output_dir {params.out_dir} \
                   -o {params.prefix}  \
                   2>> {log} 1>&2;
        """
