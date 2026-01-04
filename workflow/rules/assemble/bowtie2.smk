rule assemble__bowtie2__build_run:
    """Index a megahit assembly"""
    input:
        contigs=lambda wildcards: (
            MEGAHIT / f"{wildcards.assembly_id}.fa.gz" if config["assembler"] == "megahit" else 
            METASPADES / f"{wildcards.assembly_id}.fa.gz"if config["assembler"] == "metaspades" else 
            PROVIDED / f"{wildcards.assembly_id}.fa.gz"
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
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__bowtie2__build')['nvme']}",
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

# Test the samtools output type and set some env variables
SAMTOOLS_OUTTYPE = params["assemble"]["samtools"]["out_type"].upper()
assert SAMTOOLS_OUTTYPE in {"BAM", "CRAM"}, "samtools.out_type must be BAM or CRAM"
ALIGN_EXT = "cram" if SAMTOOLS_OUTTYPE == "CRAM" else "bam"

rule assemble__bowtie2__map:
    """Map one sample to one megahit assembly with intermediate BAM"""
    input:
        mock=ASSEMBLE_INDEX / "{assembly_id}",
        forward_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        reference=lambda wc: (
            MEGAHIT / f"{wc.assembly_id}.fa.gz" if config["assembler"] == "megahit" else
            METASPADES / f"{wc.assembly_id}.fa.gz" if config["assembler"] == "metaspades" else
            PROVIDED / f"{wc.assembly_id}.fa.gz"
        ),
    output:
        file=Path(str(ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.") + ALIGN_EXT),
        tmp_bam=temp(ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.tmp.bam"),
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
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__bowtie2__map')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__bowtie2__map"))
    params:
        samtools_mem=params["assemble"]["samtools"]["mem"],
        samtools_outtype=params["assemble"]["samtools"]["out_type"].upper(),
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    shell:
        """
        echo "=== Running assemble__bowtie2__map for assembly {wildcards.assembly_id}, sample {wildcards.sample_id} and library {wildcards.library_id} ===" > {log}.{resources.attempt} 1>&2

        # Step 1: Bowtie2 mapping
        bowtie2 \
            -x {input.mock} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
        2>> {log}.{resources.attempt} | samtools view -b -o {output.tmp_bam} - 2>> {log}.{resources.attempt}

        # Step 2: Sort and convert
        samtools sort \
            -l 9 \
            -m {params.samtools_mem} \
            -O {params.samtools_outtype} \
            -o {output.file} \
            --threads {threads} \
            {output.tmp_bam} 2>> {log}.{resources.attempt}
        
        echo "=== Finished assemble__bowtie2__map for assembly {wildcards.assembly_id}, sample {wildcards.sample_id} and library {wildcards.library_id} ===" > {log}.{resources.attempt} 1>&2

        mv {log}.{resources.attempt} {log}
        """



rule assemble__bowtie2:
    """Map all samples to all the assemblies that they belong to"""
    input:
        [
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.cram"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
