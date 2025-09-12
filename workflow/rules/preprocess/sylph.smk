rule preprocess__sylph_profile:
    """Run Sylph profiler"""
    input:
        forwards=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverses=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        db=features["databases"]["sylph"],
    output:
        PRE_SYLPH / "{sample_id}.{library_id}.profiling.tsv",
    log:
        PRE_SYLPH / "{sample_id}.{library_id}.profiling_report.log",
    conda:
        "__environment__.yml"
    container:
        docker["sylph"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"]//config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
        nvme = config["resources"]["nvme"]["small"],
    shell:
        """
        sylph profile {input.db} -1 {input.forwards} -2 {input.reverses} -t {threads} > {output} 2>> {log} 1>&2
        """


rule preprocess__sylph:
    """Run Sylph"""
    input:
        [
            PRE_SYLPH / f"{sample_id}.{library_id}.profiling.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
