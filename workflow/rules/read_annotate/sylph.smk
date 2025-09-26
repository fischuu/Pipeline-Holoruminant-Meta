rule read_annotate__sylph_profile:
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
    threads: esc("cpus", "read_annotate__sylph_profile")
    resources:
        runtime=esc("runtime", "read_annotate__sylph_profile"),
        mem_mb=esc("mem_mb", "read_annotate__sylph_profile"),
        cpus_per_task=esc("cpus", "read_annotate__sylph_profile"),
        slurm_partition=esc("partition", "read_annotate__sylph_profile"),
        #gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'read_annotate__sylph_profile')['nvme']}",
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'read_annotate__sylph_profile')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("read_annotate__sylph_profile"))
    shell:
        """
        sylph profile {input.db} -1 {input.forwards} -2 {input.reverses} -t {threads} > {output} 2> {log}
        """


rule read_annotate__sylph:
    """Run Sylph"""
    input:
        [
            PRE_SYLPH / f"{sample_id}.{library_id}.profiling.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
