rule report__preprocess:
    """
    Create the pipeline report for the preproces module (R).
    """
    input:
        rules.preprocess.input,
    output:
        html=PIPELINE_REPORT / "preprocess.html",
    log:
        PIPELINE_REPORT / "preprocess.log",
    benchmark:
        PIPELINE_REPORT / "preprocess_benchmark.tsv",
    container:
        docker["r_report"]
    params:
       script=PREPROCESS_R,
       features=config["features-file"],
       project=WD
    threads: esc("cpus", "report__preprocess")
    resources:
        runtime=esc("runtime", "report__preprocess"),
        mem_mb=esc("mem_mb", "report__preprocess"),
        cpus_per_task=esc("cpus", "report__preprocess"),
        slurm_partition=esc("partition", "report__preprocess"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'report__preprocess')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__preprocess"))
    shell:"""
       R -e "features_file <- '{params.features}'; \
             project_folder <- '{params.project}' ; \
             snakemake <- TRUE ; \
             rmarkdown::render('{params.script}',output_file=file.path('{params.project}','{output}'))" &> {log}
    """
