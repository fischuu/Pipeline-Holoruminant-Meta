rule annotate__quast:
    """Run quast over one the dereplicated mags"""
    input:
        DREP / "dereplicated_genomes.fa.gz",
    output:
        directory(QUAST),
    log:
        QUAST / "quast.log",
    container:
        docker["annotate"]
    threads: esc("cpus", "annotate__quast")
    resources:
        runtime=esc("runtime", "annotate__quast"),
        mem_mb=esc("mem_mb", "annotate__quast"),
        cpu_per_task=esc("cpus", "annotate__quast"),
        slurm_partition=esc("partition", "annotate__quast"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "annotate__quast", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("annotate__quast"))
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """
