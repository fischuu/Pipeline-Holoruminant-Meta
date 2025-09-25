rule annotate__bakta:
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
    threads: esc("cpus", "annotate__bakta")
    resources:
        runtime=esc("runtime", "annotate__bakta"),
        mem_mb=esc("mem_mb", "annotate__bakta"),
        cpus_per_task=esc("cpus", "annotate__bakta"),
        slurm_partition=esc("partition", "annotate__bakta"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'annotate__bakta')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("annotate__bakta"))
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

rule annotate__bakta_mags_run:
    """Run Bakta over the dereplicated mags"""
    input:
        contigs=MAGSCOT / "{assembly_id}.fa.gz",
    output:
        tsv=BAKTAMAG / "bakta_{assembly_id}.tsv",
        faa=BAKTAMAG / "bakta_{assembly_id}.faa",
    log:
        BAKTAMAG / "bakta_{assembly_id}.log",
    container:
        docker["annotate"]
    params:
        out_dir=BAKTAMAG,
        db=features["databases"]["bakta"],
        options=params["annotate"]["bakta"]["additional_options"],
        tmpdir=config["tmp_storage"]
    threads: esc("cpus", "annotate__bakta_mags_run")
    resources:
        runtime=esc("runtime", "annotate__bakta_mags_run"),
        mem_mb=esc("mem_mb", "annotate__bakta_mags_run"),
        cpus_per_task=esc("cpus", "annotate__bakta_mags_run"),
        slurm_partition=esc("partition", "annotate__bakta_mags_run"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'annotate__bakta_mags_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("annotate__bakta_mags_run"))
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

rule annotate__bakta_mags:
    """Run Bakta over the dereplicated mags"""
     input:
        expand(BAKTAMAG / "bakta_{assembly_id}.faa", assembly_id=ASSEMBLIES),

