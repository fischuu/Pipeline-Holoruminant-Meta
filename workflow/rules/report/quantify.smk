rule report__quantify:
    """
    Create the pipeline report for the quantify module (R).
    """
    input:
        rules.quantify.input,
    output:
        html=PIPELINE_REPORT / "quantify.html",
    log:
        PIPELINE_REPORT / "quantify.log",
    benchmark:
        PIPELINE_REPORT / "quantify_benchmark.tsv",
    container:
        docker["r_report"]
    params:
       script=QUANTIFY_R,
       features=config["features-file"],
       wd=WD
    threads: esc("cpus", "report__quantify")
    resources:
        runtime=esc("runtime", "report__quantify"),
        mem_mb=esc("mem_mb", "report__quantify"),
        cpus_per_task=esc("cpus", "report__quantify"),
        slurm_partition=esc("partition", "report__quantify"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'report__quantify')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__quantify"))
    shell:"""
       R -e "working_dir <- '{params.wd}'; \
             features_file <- '{params.features}'; \
             rmarkdown::render('{params.script}',output_file='{params.wd}/{output}')" &> {log}
    """
