rule reference__hosts__recompress:
    """Recompress the reference with bgzip"""
    input:
        fa_gz=lambda wildcards: features["hosts"][wildcards.genome],
    output:
        HOSTS / "{genome}.fa.gz",
    log:
        HOSTS / "{genome}.log",
    container:
        docker["reference"]
    threads: esc("cpus", "reference__hosts__recompress")
    resources:
        runtime=esc("runtime", "reference__hosts__recompress"),
        mem_mb=esc("mem_mb", "reference__hosts__recompress"),
        cpus_per_task=esc("cpus", "reference__hosts__recompress"),
        slurm_partition=esc("partition", "reference__hosts__recompress"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "reference__hosts__recompress", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__hosts__recompress"))
    shell:
        """
        echo "$(date) **Starting rule reference__hosts__recompress**" > {log}

        (gzip -dc {input.fa_gz} | bgzip -@ {threads} > {output}) 2> {log}
        
        echo "$(date) **Finished rule reference__hosts__recompress**" >> {log}
        """


rule reference__hosts:
    """Run all the steps for reference(s) preparations"""
    input:
        [HOSTS / f"{genome}.fa.gz" for genome in HOST_NAMES],
