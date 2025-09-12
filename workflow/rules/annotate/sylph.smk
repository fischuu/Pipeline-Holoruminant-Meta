rule _annotate__sylph_profile:
    """Run Sylph profiler"""
    input:
        mags=DREP / "dereplicated_genomes",
        db=features["databases"]["sylph"],
    output:
        SYLPH / "profiling.tsv",
    log:
        SYLPH / "quality_report.log",
    container:
        docker["sylph"]
    params:
        out_dir=CHECKM / "predict",
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"]//config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
        nvme = config["resources"]["nvme"]["small"],
    shell:
        """
        sylph profile {input.db} {input.mags}/* -t {threads} > {output} 2>> {log} 1>&2
        """


rule annotate__sylph:
    """Run Sylph"""
    input:
        rules._annotate__sylph_profile.output,
