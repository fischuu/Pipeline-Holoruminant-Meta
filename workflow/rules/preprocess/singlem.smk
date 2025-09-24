rule preprocess__singlem__pipe:
    """Run singlem over one sample

    Note: SingleM asks in the documentation for the raw reads. Here we are
    passing it the non-host and trimmed ones.
    """
    input:
        forward_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        data=features["databases"]["singlem"], 
    output:
        archive_otu_table=SINGLEM / "pipe" / "{sample_id}.{library_id}.archive.json",
        otu_table=SINGLEM / "pipe" / "{sample_id}.{library_id}.otu_table.tsv",
        condense=SINGLEM / "pipe" / "{sample_id}.{library_id}.condense.tsv",
    log:
        SINGLEM / "pipe" / "{sample_id}.{library_id}.log",
    benchmark:
        SINGLEM / "benchmark/pipe" / "{sample_id}.{library_id}.tsv",
    container:
        docker["preprocess"]
    threads: esc("cpus", "preprocess__singlem__pipe")
    resources:
        runtime=esc("runtime", "preprocess__singlem__pipe"),
        mem_mb=esc("mem_mb", "preprocess__singlem__pipe"),
        cpu_per_task=esc("cpus", "preprocess__singlem__pipe"),
        slurm_partition=esc("partition", "preprocess__singlem__pipe"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__singlem__pipe", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__singlem__pipe"))
    params:
        tmp = config["tmp_storage"]
    shell:
        """
        TMPDIR={params.tmp}
        
        echo "Checking disk space for TMPDIR: ${{TMPDIR:-/tmp}}" >> {log}
        df -h ${{TMPDIR:-/tmp}} >> {log}
        echo "Disk space check completed." >> {log}

        singlem pipe \
            --forward {input.forward_} \
            --reverse {input.reverse_} \
            --otu-table {output.otu_table} \
            --archive-otu-table {output.archive_otu_table} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.data} \
            --threads {threads} \
            --assignment-threads {threads} \
        2>> {log} 1>&2 || true
        """


rule preprocess__singlem__condense:
    """Aggregate all the singlem results into a single table"""
    input:
        archive_otu_tables=[
            SINGLEM / "pipe" / f"{sample_id}.{library_id}.archive.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        database=features["databases"]["singlem"],
    output:
        condense=SINGLEM / "singlem.tsv",
    log:
        SINGLEM / "singlem.log",
    benchmark:
        SINGLEM / "benchmark/singlem.tsv",
    container:
        docker["preprocess"]
    params:
        input_dir=SINGLEM,
    threads: esc("cpus", "preprocess__singlem__condense")
    resources:
        runtime=esc("runtime", "preprocess__singlem__condense"),
        mem_mb=esc("mem_mb", "preprocess__singlem__condense"),
        cpu_per_task=esc("cpus", "preprocess__singlem__condense"),
        slurm_partition=esc("partition", "preprocess__singlem__condense"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__singlem__condense", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__singlem__condense"))
    shell:
        """
        singlem condense \
            --input-archive-otu-tables {input.archive_otu_tables} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.database} \
        2> {log} 1>&2
        """


rule preprocess__singlem__microbial_fraction:
    """Run singlem microbial_fraction over one sample"""
    input:
        forward_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        data=features["databases"]["singlem"],
        condense=SINGLEM / "pipe" / "{sample_id}.{library_id}.condense.tsv",
    output:
        microbial_fraction=SINGLEM
        / "microbial_fraction"
        / "{sample_id}.{library_id}.tsv",
    log:
        SINGLEM / "microbial_fraction" / "{sample_id}.{library_id}.log",
    benchmark:
        SINGLEM / "benchmark/microbial_fraction" / "{sample_id}.{library_id}.tsv"
    container:
        docker["preprocess"]
    threads: esc("cpus", "preprocess__singlem__microbial_fraction")
    resources:
        runtime=esc("runtime", "preprocess__singlem__microbial_fraction"),
        mem_mb=esc("mem_mb", "preprocess__singlem__microbial_fraction"),
        cpu_per_task=esc("cpus", "preprocess__singlem__microbial_fraction"),
        slurm_partition=esc("partition", "preprocess__singlem__microbial_fraction"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__singlem__microbial_fraction", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__singlem__microbial_fraction"))
    shell:
        """
        singlem microbial_fraction \
            --forward {input.forward_} \
            --reverse {input.reverse_} \
            --input-profile {input.condense} \
            --output-tsv {output.microbial_fraction} \
            --metapackage {input.data} \
        2> {log} 1>&2
        """


rule preprocess__singlem__aggregate_microbial_fraction:
    """Aggregate all the microbial_fraction files into one tsv"""
    input:
        tsvs=[
            SINGLEM / "microbial_fraction" / f"{sample_id}.{library_id}.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        tsv=SINGLEM / "microbial_fraction.tsv",
    log:
        SINGLEM / "microbial_fraction.log",
    benchmark:
        SINGLEM / "benchmark/microbial_fraction.tsv"
    container:
        docker["preprocess"]
    threads: esc("cpus", "preprocess__singlem__aggregate_microbial_fraction")
    resources:
        runtime=esc("runtime", "preprocess__singlem__aggregate_microbial_fraction"),
        mem_mb=esc("mem_mb", "preprocess__singlem__aggregate_microbial_fraction"),
        cpu_per_task=esc("cpus", "preprocess__singlem__aggregate_microbial_fraction"),
        slurm_partition=esc("partition", "preprocess__singlem__aggregate_microbial_fraction"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__singlem__aggregate_microbial_fraction", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__singlem__aggregate_microbial_fraction"))
    shell:
        """
        ( csvstack \
            --tabs \
            {input.tsvs} \
        | csvformat \
            --out-tabs \
        > {output.tsv} \
        ) 2> {log}
        """

rule preprocess__singlem:
    input:
        rules.preprocess__singlem__condense.output,
        rules.preprocess__singlem__aggregate_microbial_fraction.output,
