rule assemble__bowtie2__build_run:
    """Index a megahit assembly"""
    input:
        contigs=lambda wildcards: (
            MEGAHIT / f"{wildcards.assembly_id}.fa.gz" if config["assembler"] == "megahit" else
            METASPADES / f"{wildcards.assembly_id}.fa.gz"
        )
    output:
        mock=touch(ASSEMBLE_INDEX / "{assembly_id}"),
    log:
        ASSEMBLE_INDEX / "{assembly_id}.log",
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__bowtie2__build")
    resources:
        runtime=esc("runtime", "assemble__bowtie2__build"),
        mem_mb=esc("mem_mb", "assemble__bowtie2__build"),
        cpus_per_task=esc("cpus", "assemble__bowtie2__build"),
        slurm_partition=esc("partition", "assemble__bowtie2__build"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "assemble__bowtie2__build", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__bowtie2__build"))
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.contigs} \
            {output.mock} \
        2> {log}.{resources.attempt} 1>&2

        mv {log}.{resources.attempt} {log}
        """


rule assemble__bowtie2__build:
    """Index all megahit assemblies"""
    input:
        [ASSEMBLE_INDEX / f"{assembly_id}" for assembly_id in ASSEMBLIES],


rule assemble__bowtie2__map:
    """Map one sample to one megahit assembly"""
    input:
        mock=ASSEMBLE_INDEX / "{assembly_id}",
        forward_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        reference=lambda wildcards: (
            MEGAHIT / f"{wildcards.assembly_id}.fa.gz" if config["assembler"] == "megahit" else
            METASPADES / f"{wildcards.assembly_id}.fa.gz"
        ),
        fai=lambda wildcards: (
            MEGAHIT / f"{wildcards.assembly_id}.fa.gz.fai" if config["assembler"] == "megahit" else
            METASPADES / f"{wildcards.assembly_id}.fa.gz.fai"
        ),
    output:
        cram=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram",
    log:
        log=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.log",
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__bowtie2__map")
    resources:
        runtime=esc("runtime", "assemble__bowtie2__map"),
        mem_mb=esc("mem_mb", "assemble__bowtie2__map"),
        cpus_per_task=esc("cpus", "assemble__bowtie2__map"),
        slurm_partition=esc("partition", "assemble__bowtie2__map"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "assemble__bowtie2__map", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__bowtie2__map"))
    params:
        samtools_mem=params["assemble"]["samtools"]["mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
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
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2>> {log}.{resources.attempt} 1>&2

        mv {log}.{resources.attempt} {log}
        """


rule assemble__bowtie2:
    """Map all samples to all the assemblies that they belong to"""
    input:
        [
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.cram"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
