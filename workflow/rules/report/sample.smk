rule _report__sample__multiqc:
    input:
        get_stats_files_from_sample_and_library_ids,
    output:
        html=REPORT_SAMPLE / "{sample_id}.{library_id}.html",
    log:
        REPORT_SAMPLE / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["report"]
    params:
        dir=REPORT_SAMPLE,
        filename=lambda wildcards: f"{wildcards.sample_id}.{wildcards.library_id}",
    resources:
        mem_mb=8 * 1024,
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
