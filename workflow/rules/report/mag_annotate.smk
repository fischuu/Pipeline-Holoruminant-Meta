rule report__mag_annotate:
    """
    Create the pipeline report for the mag_annotate module (R).
    """
    input:
        rules.mag_annotate.input,
    output:
        html=PIPELINE_REPORT / "mag_annotate.html",
    log:
        PIPELINE_REPORT / "mag_annotate.log",
    benchmark:
        PIPELINE_REPORT / "mag_annotate_benchmark.tsv",
    container:
        docker["r_report"]
    params:
       script=MAG_ANNOTATE_R,
       features=config["features-file"],
       wd=WD
    threads: esc("cpus", "report__mag_annotate")
    resources:
        runtime=esc("runtime", "report__mag_annotate"),
        mem_mb=esc("mem_mb", "report__mag_annotate"),
        cpus_per_task=esc("cpus", "report__mag_annotate"),
        slurm_partition=esc("partition", "report__mag_annotate"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'report__mag_annotate')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__mag_annotate"))
    shell:"""
       R -e "working_dir <- '{params.wd}'; \
             features_file <- '{params.features}'; \
             rmarkdown::render('{params.script}',output_file='{params.wd}/{output}')" &> {log}
    """
