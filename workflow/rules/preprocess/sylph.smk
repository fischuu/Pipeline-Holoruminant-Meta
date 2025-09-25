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
    container:
        docker["sylph"]
    threads: esc("cpus", "preprocess__sylph_profile")
    resources:
        runtime=esc("runtime", "preprocess__sylph_profile"),
        mem_mb=esc("mem_mb", "preprocess__sylph_profile"),
        cpus_per_task=esc("cpus", "preprocess__sylph_profile"),
        slurm_partition=esc("partition", "preprocess__sylph_profile"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__sylph_profile", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__sylph_profile"))
    shell:
        """
        sylph profile {input.db} -1 {input.forwards} -2 {input.reverses} -t {threads} > {output} 2> {log}
        """


rule preprocess__sylph:
    """Run Sylph"""
    input:
        [
            PRE_SYLPH / f"{sample_id}.{library_id}.profiling.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
