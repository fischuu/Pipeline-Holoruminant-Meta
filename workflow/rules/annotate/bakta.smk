rule _annotate__bakta:
    """Run Bakta over the dereplicated mags"""
    input:
        contigs=DREP / "dereplicated_genomes.fa.gz",
    output:
        tsv=BAKTA / "bakta.tsv",
        faa=BAKTA / "bakta.faa",
    log:
        BAKTA / "bakta.log",
    conda:
        "annotate.yml"
    singularity:
        docker["annotate"]
    params:
        out_dir=BAKTA,
        db=features["databases"]["bakta"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"]//config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        bakta --db {params.db} \
              --force \
              --output {params.out_dir} \
              --prefix bakta \
              --threads {threads} \
              {input.contigs} \
              2>> {log} 1>&2
        """
