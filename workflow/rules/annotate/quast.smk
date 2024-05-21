rule annotate__quast:
    """Run quast over one the dereplicated mags"""
    input:
        DREP / "dereplicated_genomes.fa.gz",
    output:
        directory(QUAST),
    log:
        QUAST / "quast.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["annotate"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """
