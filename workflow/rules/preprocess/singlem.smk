rule _preprocess__singlem__pipe:
    """Run singlem over one sample

    Note: SingleM asks in the documentation for the raw reads. Here we are
    passing it the non-host and trimmed ones.
    """
    input:
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        data=features["databases"]["singlem"],
    output:
        archive_otu_table=SINGLEM / "pipe" / "{sample_id}.{library_id}.archive.json",
        otu_table=SINGLEM / "pipe" / "{sample_id}.{library_id}.otu_table.tsv",
        condense=SINGLEM / "pipe" / "{sample_id}.{library_id}.condense.tsv",
    log:
        SINGLEM / "pipe" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
    threads: 1
    resources:
        runtime=48 * 60,
        mem_mb=16 * 1024,
    shell:
        """
        singlem pipe \
            --forward {input.forward_} \
            --reverse {input.reverse_} \
            --otu-table {output.otu_table} \
            --archive-otu-table {output.archive_otu_table} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.data} \
            --threads {threads} \
            --assignment-threads {threads} \
        2> {log} 1>&2 || true
        """


rule _preprocess__singlem__condense:
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
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
    params:
        input_dir=SINGLEM,
    resources:
        runtime=6 * 60,
    shell:
        """
        singlem condense \
            --input-archive-otu-tables {input.archive_otu_tables} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.database} \
        2> {log} 1>&2
        """


rule _preprocess__singlem__microbial_fraction:
    """Run singlem microbial_fraction over one sample"""
    input:
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        data=features["databases"]["singlem"],
        condense=SINGLEM / "pipe" / "{sample_id}.{library_id}.condense.tsv",
    output:
        microbial_fraction=SINGLEM
        / "microbial_fraction"
        / "{sample_id}.{library_id}.tsv",
    log:
        SINGLEM / "microbial_fraction" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
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


rule _preprocess__singlem__aggregate_microbial_fraction:
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
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
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
        rules._preprocess__singlem__condense.output,
        rules._preprocess__singlem__aggregate_microbial_fraction.output,
