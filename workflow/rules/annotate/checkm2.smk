rule annotate__checkm2__predict:
    """Run CheckM2 over the dereplicated mags"""
    input:
        mags=DREP / "dereplicated_genomes",
        db=features["databases"]["checkm2"],
    output:
        CHECKM / "quality_report.tsv",
    log:
        CHECKM / "quality_report.log",
    container:
        docker["checkm2"]
    params:
        out_dir=CHECKM / "predict",
    threads: esc("cpus", "annotate__checkm2__predict")
    resources:
        runtime=esc("runtime", "annotate__checkm2__predict"),
        mem_mb=esc("mem_mb", "annotate__checkm2__predict"),
        cpu_per_task=esc("cpus", "annotate__checkm2__predict"),
        slurm_partition=esc("partition", "annotate__checkm2__predict"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "annotate__checkm2__predict", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("annotate__checkm2__predict"))
    shell:
        """
        rm -rfv {params.out_dir} 2> {log} 1>&2

        checkm2 predict \
            --threads {threads} \
            --input {input.mags} \
            --extension .fa.gz \
            --output-directory {params.out_dir} \
            --database_path {input.db} \
            --remove_intermediates \
        2>> {log} 1>&2

        mv {params.out_dir}/quality_report.tsv {output} 2>> {log} 1>&2
        rm --recursive --verbose --force {params.out_dir} 2>> {log} 1>&2
        """


rule annotate__checkm2:
    """Run CheckM2"""
    input:
        rules.annotate__checkm2__predict.output,
