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
    conda:
        "__environment__.yml"
    container:
        docker["r_report"]
    params:
       script=READS_R,
       features=config["features-file"],
       wd=WD
    threads: config["resources"]["cpu_per_task"]["single_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["single_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["lowmem"],
        time=config["resources"]["time"]["shortrun"],
    shell:"""
       R -e "working_dir <- '{params.wd}'; \
             features_file <- '{params.features}'; \
             rmarkdown::render('{params.script}',output_file='{params.wd}/{output}')" &> {log}
    """