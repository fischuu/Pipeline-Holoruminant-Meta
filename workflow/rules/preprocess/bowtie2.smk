rule _preprocess__bowtie2__build:
    """Build PRE_BOWTIE2 index for the host reference

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=HOSTS / "{genome}.fa.gz",
        faidx=HOSTS / "{genome}.fa.gz.fai",
    output:
        mock=touch(PRE_INDEX / "{genome}"),
    log:
        PRE_INDEX / "{genome}.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
        attempt=get_attempt,
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.reference} \
            {output.mock} \
        2> {log}.{resources.attempt} 1>&2

        mv \
            {log}.{resources.attempt} \
            {log}
        """


rule _preprocess__bowtie2__map:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=get_input_forward_for_host_mapping,
        reverse_=get_input_reverse_for_host_mapping,
        mock=PRE_INDEX / "{genome}",
        reference=HOSTS / "{genome}.fa.gz",
        faidx=HOSTS / "{genome}.fa.gz.fai",
    output:
        cram=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram",
    log:
        PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.log",
    params:
        samtools_mem=params["preprocess"]["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
        attempt=get_attempt,
    shell:
        """
        find \
            $(dirname {output.cram}) \
            -name "$(basename {output.cram}).tmp.*.bam" \
            -delete \
        2> {log}.{resources.attempt} 1>&2

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
        ) 2>> {log}.{resources.attempt} 1>&2

        mv \
            {log}.{resources.attempt} \
            {log}
        """


rule _preprocess__bowtie2__extract_nonhost:
    """
    Keep only pairs unmapped to the human reference genome, sort by name rather
    than by coordinate, and convert to FASTQ.
    """
    input:
        cram=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram",
        reference=HOSTS / "{genome}.fa.gz",
        fai=HOSTS / "{genome}.fa.gz.fai",
    output:
        forward_=PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}_2.fq.gz",
    log:
        PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
    params:
        samtools_mem=params["preprocess"]["bowtie2"]["samtools"]["mem_per_thread"],
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"]
    shell:
        """
        ( samtools view \
            --reference {input.reference} \
            --threads {threads} \
            -u \
            -o /dev/stdout \
            -f 12 \
            {input.cram} \
        | samtools collate \
            -O \
            -u \
            -f \
            --reference {input.reference} \
            -@ {threads} \
            - \
        | samtools fastq \
            -1 >(pigz > {output.forward_}) \
            -2 >(pigz > {output.reverse_}) \
            -0 /dev/null \
            -c 9 \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule preprocess__bowtie2__extract_nonhost:
    """Run bowtie2_extract_nonchicken_one for all PE libraries"""
    input:
        [
            PRE_BOWTIE2 / f"non{genome}" / f"{sample_id}.{library_id}_{end}.fq.gz"
            for genome in [LAST_HOST]
            if LAST_HOST
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule preprocess__bowtie2:
    """Run all the preprocessing steps for bowtie2"""
    input:
        rules.preprocess__bowtie2__extract_nonhost.input,
