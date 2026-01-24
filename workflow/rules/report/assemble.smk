rule report__assemble:
    """
    Create the pipeline report for the assemble module (R).
    """
    input:
        rules.assemble.input,
    output:
        html=PIPELINE_REPORT / "assemble.html",
    log:
        PIPELINE_REPORT / "assemble.log",
    benchmark:
        PIPELINE_REPORT / "assemble_benchmark.tsv",
    container:
        docker["r_report"]
    params:
       script=ASSEMBLE_R,
       features=config["features-file"],
       wd=WD
    threads: esc("cpus", "report__assemble")
    resources:
        runtime=esc("runtime", "report__assemble"),
        mem_mb=esc("mem_mb", "report__assemble"),
        cpus_per_task=esc("cpus", "report__assemble"),
        slurm_partition=esc("partition", "report__assemble"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'report__assemble')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__assemble"))
    shell:"""
       R -e "working_dir <- '{params.wd}'; \
             features_file <- '{params.features}'; \
             rmarkdown::render('{params.script}',output_file='{params.wd}/{output}')" &> {log}
    """
