rule report__reads:
    """
    Create the pipeline report for the reads module (R).
    """
    input:
        rules.reads.input,
    output:
        html=PIPELINE_REPORT / "reads.html",
    log:
        PIPELINE_REPORT / "reads.log",
    benchmark:
        PIPELINE_REPORT / "reads_benchmark.tsv",
    container:
        docker["r_report"]
    params:
       script=READS_R,
       features=config["features-file"],
       wd=WD
    threads: esc("cpus", "report__reads")
    resources:
        runtime=esc("runtime", "report__reads"),
        mem_mb=esc("mem_mb", "report__reads"),
        cpus_per_task=esc("cpus", "report__reads"),
        slurm_partition=esc("partition", "report__reads"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'report__reads')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__reads"))
    shell:"""
       R -e "working_dir <- '{params.wd}'; \
             features_file <- '{params.features}'; \
             rmarkdown::render('{params.script}',output_file='{params.wd}/{output}')" &> {log}
    """
