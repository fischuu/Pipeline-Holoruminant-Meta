rule report__read_annotate:
    """
    Create the pipeline report for the read_annotate module (R).
    """
    input:
        rules.read_annotate.input,
    output:
        html=PIPELINE_REPORT / "read_annotate.html",
    log:
        PIPELINE_REPORT / "read_annotate.log",
    benchmark:
        PIPELINE_REPORT / "read_annotate_benchmark.tsv",
    container:
        docker["r_report"]
    params:
       script=READ_ANNOTATE_R,
       features=config["features-file"],
       wd=WD
    threads: esc("cpus", "report__read_annotate")
    resources:
        runtime=esc("runtime", "report__read_annotate"),
        mem_mb=esc("mem_mb", "report__read_annotate"),
        cpus_per_task=esc("cpus", "report__read_annotate"),
        slurm_partition=esc("partition", "report__read_annotate"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'report__read_annotate')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__read_annotate"))
    shell:"""
       R -e "working_dir <- '{params.wd}'; \
             features_file <- '{params.features}'; \
             rmarkdown::render('{params.script}',output_file='{params.wd}/{output}')" &> {log}
    """
