rule report__sample__multiqc:
    input:
        get_stats_files_from_sample_and_library_ids,
    output:
        html=REPORT_SAMPLE / "{sample_id}.{library_id}.html",
    log:
        REPORT_SAMPLE / "{sample_id}.{library_id}.log",
    container:
        docker["report"]
    params:
        dir=REPORT_SAMPLE,
        filename=lambda wildcards: f"{wildcards.sample_id}.{wildcards.library_id}",
    threads: esc("cpus", "report__sample__multiqc")
    resources:
        runtime=esc("runtime", "report__sample__multiqc"),
        mem_mb=esc("mem_mb", "report__sample__multiqc"),
        cpu_per_task=esc("cpus", "report__sample__multiqc"),
        slurm_partition=esc("partition", "report__sample__multiqc"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "report__sample__multiqc", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__sample__multiqc"))
    shell:
        """
        multiqc \
            --filename {params.filename} \
            --title {params.filename} \
            --force \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__sample:
    input:
        [
            REPORT_SAMPLE / f"{sample_id}.{library_id}.html"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
