rule report__contig_annotate:
    """
    Create the pipeline report for the contig_annotate module (R).
    """
    input:
        rules.contig_annotate.input,
    output:
        html=PIPELINE_REPORT / "contig_annotate.html",
    log:
        PIPELINE_REPORT / "contig_annotate.log",
    benchmark:
        PIPELINE_REPORT / "contig_annotate_benchmark.tsv",
    container:
        docker["r_report"]
    params:
       script=CONTIG_ANNOTATE_R,
       features=config["features-file"],
       wd=WD
    threads: esc("cpus", "report__contig_annotate")
    resources:
        runtime=esc("runtime", "report__contig_annotate"),
        mem_mb=esc("mem_mb", "report__contig_annotate"),
        cpus_per_task=esc("cpus", "report__contig_annotate"),
        slurm_partition=esc("partition", "report__contig_annotate"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'report__contig_annotate')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__contig_annotate"))
    shell:"""
       R -e "working_dir <- '{params.wd}'; \
             features_file <- '{params.features}'; \
             rmarkdown::render('{params.script}',output_file='{params.wd}/{output}')" &> {log}
    """
