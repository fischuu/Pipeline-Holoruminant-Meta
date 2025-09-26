rule read_annotate__nonpareil__run:
    """Run nonpareil over one sample

    Note: Nonpareil only ask for one of the pair-end reads
    Note2: it has to be fastq. The process substitution trick does not work
    Note3: in case that nonpareil fails for low coverage samples, it creates
    empty files
    """
    input:
        forward_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
    output:
        npa=touch(NONPAREIL / "{sample_id}.{library_id}.npa"),
        npc=touch(NONPAREIL / "{sample_id}.{library_id}.npc"),
        npl=touch(NONPAREIL / "{sample_id}.{library_id}.npl"),
        npo=touch(NONPAREIL / "{sample_id}.{library_id}.npo"),
    log:
        NONPAREIL / "{sample_id}.{library_id}.log",
    benchmark:
        NONPAREIL / "benchmark/{sample_id}.{library_id}.tsv",
    container:
        docker["preprocess"]
    params:
        prefix=lambda w: NONPAREIL / f"{w.sample_id}.{w.library_id}",
        reads=lambda w: NONPAREIL /  f"{w.sample_id}.{w.library_id}_1.fq",
        X=params["read_annotate"]["nonpareil"]["X"],
        tmp = config["tmp_storage"]
    threads: esc("cpus", "read_annotate__nonpareil__run")
    resources:
        runtime=esc("runtime", "read_annotate__nonpareil__run"),
        mem_mb=esc("mem_mb", "read_annotate__nonpareil__run"),
        cpus_per_task=esc("cpus", "read_annotate__nonpareil__run"),
        slurm_partition=esc("partition", "read_annotate__nonpareil__run"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'read_annotate__nonpareil__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("read_annotate__nonpareil__run"))
    shell:
        """
        TMPDIR={params.tmp}
        
        gunzip -f -c {input.forward_} > {params.reads} 2> {log}

        nonpareil \
            -s {params.reads} \
            -T kmer \
            -b {params.prefix} \
            -f fastq \
            -t {threads} \
            -X {params.X} \
        2>> {log} 1>&2

        rm --force {params.reads} 2>> {log} 1>&2
        """


rule read_annotate__nonpareil__aggregate:
    """Aggregate all the nonpareil results into a single table"""
    input:
        [
            NONPAREIL / f"{sample_id}.{library_id}.npo"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        NONPAREIL / "nonpareil.tsv",
    log:
        NONPAREIL / "nonpareil.log",
    benchmark:
        NONPAREIL / "benchmark/nonpareil.tsv",
    container:
        docker["preprocess"]
    threads: esc("cpus", "read_annotate__nonpareil__aggregate")
    resources:
        runtime=esc("runtime", "read_annotate__nonpareil__aggregate"),
        mem_mb=esc("mem_mb", "read_annotate__nonpareil__aggregate"),
        cpus_per_task=esc("cpus", "read_annotate__nonpareil__aggregate"),
        slurm_partition=esc("partition", "read_annotate__nonpareil__aggregate"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'read_annotate__nonpareil__aggregate')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("read_annotate__nonpareil__aggregate"))
    params:
        input_dir=NONPAREIL,
        script_folder=SCRIPT_FOLDER,
    shell:
        """
        Rscript --no-init-file {params.script_folder}/aggregate_nonpareil.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule read_annotate__nonpareil:
    input:
        rules.read_annotate__nonpareil__aggregate.output,
