rule _report__preprocess:
    """
    Create the pipeline report for the preproces module (R).
    """
    input:
        rules.preprocess.input,
    output:
        html=REPORT_STEP / "preprocess.html",
    log:
        REPORT_STEP / "preprocess.log",
    benchmark:
        REPORT_STEP / "preprocess_benchmark.tsv",
    conda:
        "__environment__.yml"
    singularity:
        docker["report"]
    params:
       script=PREPROCESS_R,
       features=config["features-file"],
    threads: config["resources"]["cpu_per_task"]["single_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["single_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time=config["resources"]["time"]["longrun"],
    shell:"""
       R -e "features_file <- '{params.features}'; \
             rmarkdown::render('{params.script}',output_file='{output}')" &> {log}
    """