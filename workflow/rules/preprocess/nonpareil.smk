rule _preprocess__nonpareil__run:
    """Run nonpareil over one sample

    Note: Nonpareil only ask for one of the pair-end reads
    Note2: it has to be fastq. The process substitution trick does not work
    Note3: in case that nonpareil fails for low coverage samples, it creates
    empty files
    """
    input:
        forward_=get_final_forward_from_pre,
    output:
        npa=touch(NONPAREIL / "{sample_id}.{library_id}.npa"),
        npc=touch(NONPAREIL / "{sample_id}.{library_id}.npc"),
        npl=touch(NONPAREIL / "{sample_id}.{library_id}.npl"),
        npo=touch(NONPAREIL / "{sample_id}.{library_id}.npo"),
    log:
        NONPAREIL / "{sample_id}.{library_id}.log",
    benchmark:
        NONPAREIL / "benchmark/{sample_id}.{library_id}.tsv",
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
    params:
        prefix=lambda w: NONPAREIL / f"{w.sample_id}.{w.library_id}",
        reads=lambda w: NONPAREIL /  f"{w.sample_id}.{w.library_id}_1.fq",
        X=params["preprocess"]["nonpareil"]["X"],
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"]
    shell:
        """
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


rule _preprocess__nonpareil__aggregate:
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
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"]
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


rule preprocess__nonpareil:
    input:
        rules._preprocess__nonpareil__aggregate.output,
