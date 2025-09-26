rule mag_annotate__bakta:
    """Run Bakta over the dereplicated mags"""
    input:
        contigs=DREP / "dereplicated_genomes.fa.gz",
    output:
        tsv=BAKTA / "bakta.tsv",
        faa=BAKTA / "bakta.faa",
    log:
        BAKTA / "bakta.log",
    container:
        docker["mag_annotate"]
    params:
        out_dir=BAKTA,
        db=features["databases"]["bakta"],
        options=params["mag_annotate"]["bakta"]["additional_options"],
        tmpdir=config["tmp_storage"]
    threads: esc("cpus", "mag_annotate__bakta")
    resources:
        runtime=esc("runtime", "mag_annotate__bakta"),
        mem_mb=esc("mem_mb", "mag_annotate__bakta"),
        cpus_per_task=esc("cpus", "mag_annotate__bakta"),
        slurm_partition=esc("partition", "mag_annotate__bakta"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'mag_annotate__bakta')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("mag_annotate__bakta"))
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

rule mag_annotate__bakta_mags_run:
    """Run Bakta over the dereplicated mags"""
    input:
        contigs=MAGSCOT / "{assembly_id}.fa.gz",
    output:
        tsv=BAKTAMAG / "bakta_{assembly_id}.tsv",
        faa=BAKTAMAG / "bakta_{assembly_id}.faa",
    log:
        BAKTAMAG / "bakta_{assembly_id}.log",
    container:
        docker["mag_annotate"]
    params:
        out_dir=BAKTAMAG,
        db=features["databases"]["bakta"],
        options=params["mag_annotate"]["bakta"]["additional_options"],
        tmpdir=config["tmp_storage"]
    threads: esc("cpus", "mag_annotate__bakta_mags_run")
    resources:
        runtime=esc("runtime", "mag_annotate__bakta_mags_run"),
        mem_mb=esc("mem_mb", "mag_annotate__bakta_mags_run"),
        cpus_per_task=esc("cpus", "mag_annotate__bakta_mags_run"),
        slurm_partition=esc("partition", "mag_annotate__bakta_mags_run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'mag_annotate__bakta_mags_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("mag_annotate__bakta_mags_run"))
    shell:
        """
        bakta --db {params.db} \
              --force \
              --verbose \
              --output {params.out_dir} \
              --prefix bakta_{wildcards.assembly_id} \
              --threads {threads} \
              {params.options} \
              {input.contigs} \
              > {log} 2>&1
        """

rule mag_annotate__bakta_mags:
    """Run Bakta over the dereplicated mags"""
     input:
        expand(BAKTAMAG / "bakta_{assembly_id}.faa", assembly_id=ASSEMBLIES),

