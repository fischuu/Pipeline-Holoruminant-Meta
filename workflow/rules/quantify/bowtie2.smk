rule _quantify__bowtie2__build:
    """Index dereplicader"""
    input:
        contigs=DREP / "dereplicated_genomes.fa.gz",
    output:
        mock=touch(QUANT_INDEX / "dereplicated_genomes"),
    log:
        QUANT_INDEX / "dereplicated_genomes.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        mem_mb=double_ram(params["quantify"]["bowtie2-build"]["memory_gb"]),
        runtime=24 * 60,
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.contigs} \
            {output.mock} \
        2> {log} 1>&2
        """


rule _quantify__bowtie2__map:
    """Align one sample to the dereplicated genomes"""
    input:
        mock=QUANT_INDEX / "dereplicated_genomes",
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        reference=DREP / "dereplicated_genomes.fa.gz",
        fai=DREP / "dereplicated_genomes.fa.gz.fai",
    output:
        cram=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram",
    log:
        QUANT_BOWTIE2 / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    threads: 24
    params:
        samtools_mem=params["quantify"]["samtools"]["mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    resources:
        mem_mb=double_ram(params["quantify"]["bowtie2"]["memory_gb"]),
        runtime=24 * 60,
    shell:
        """
        find \
            $(dirname {output.cram}) \
            -name "$(basename {output.cram}).tmp.*.bam" \
            -delete \
        2> {log} 1>&2

        ( bowtie2 \
            -x {input.mock} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
        | samtools sort \
            -l 9 \
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2>> {log} 1>&2
        """


rule quantify__bowtie2:
    """Align all samples to the dereplicated genomes"""
    input:
        [
            QUANT_BOWTIE2 / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
