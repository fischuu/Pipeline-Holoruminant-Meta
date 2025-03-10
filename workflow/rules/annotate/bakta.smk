rule _annotate__bakta:
    """Run Bakta over the dereplicated mags"""
    input:
        contigs=DREP / "dereplicated_genomes.fa.gz",
    output:
        tsv=BAKTA / "bakta.tsv",
        faa=BAKTA / "bakta.faa",
    log:
        BAKTA / "bakta.log",
    container:
        docker["annotate"]
    params:
        out_dir=BAKTA,
        db=features["databases"]["bakta"],
        options=params["annotate"]["bakta"]["additional_options"],
        tmpdir=config["tmp_storage"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"]//config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["verylongrun"],
        nvme = config["resources"]["nvme"]["small"],
        partition = config["resources"]["partition"]["longrun"]
    shell:
        """
        bakta --db {params.db} \
              --force \
              --verbose \
              --output {params.out_dir} \
              --prefix bakta \
              --threads {threads} \
              {params.options} \
              {input.contigs} \
              > {log} 2>&1
        """
