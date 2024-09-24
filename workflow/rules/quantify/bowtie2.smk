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
    container:
        docker["quantify"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
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
    container:
        docker["quantify"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
    params:
        samtools_mem=params["quantify"]["samtools"]["mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
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
