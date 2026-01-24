rule report__reference:
    """
    Create the pipeline report for the reference module (R).
    """
    input:
        rules.reference.input,
    output:
        html=PIPELINE_REPORT / "reference.html",
    log:
        PIPELINE_REPORT / "reference.log",
    benchmark:
        PIPELINE_REPORT / "reference_benchmark.tsv",
    container:
        docker["r_report"]
    params:
       script=REFERENCE_R,
       features=config["features-file"],
       wd=WD
    threads: esc("cpus", "report__reference")
    resources:
        runtime=esc("runtime", "report__reference"),
        mem_mb=esc("mem_mb", "report__reference"),
        cpus_per_task=esc("cpus", "report__reference"),
        slurm_partition=esc("partition", "report__reference"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'report__reference')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__reference"))
    shell:"""
       R -e "working_dir <- '{params.wd}'; \
             features_file <- '{params.features}'; \
             rmarkdown::render('{params.script}',output_file='{params.wd}/{output}')" &> {log}
    """
