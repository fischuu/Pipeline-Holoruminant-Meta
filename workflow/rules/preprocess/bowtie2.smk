# ---------- preprocess__bowtie2__build ----------
rule preprocess__bowtie2__build:
    """Build PRE_BOWTIE2 index for the host reference"""
    input:
        reference=HOSTS / "{genome}.fa.gz",
        faidx=HOSTS / "{genome}.fa.gz.fai",
    output:
        mock=touch(PRE_INDEX / "{genome}"),
    log:
        PRE_INDEX / "{genome}.log",
    benchmark:
        PRE_INDEX / "benchmark/{genome}.tsv",
    container:
        docker["bowtie2"]
    threads: esc("cpus", "preprocess__bowtie2__build")
    resources:
        runtime=esc("runtime", "preprocess__bowtie2__build"),
        mem_mb=esc("mem_mb", "preprocess__bowtie2__build"),
        cpus_per_task=esc("cpus", "preprocess__bowtie2__build"),
        partition=esc("partition", "preprocess__bowtie2__build"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'preprocess__bowtie2__build')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__bowtie2__build"))
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.reference} \
            {output.mock} \
        2> {log}.{resources.attempt} 1>&2

        mv {log}.{resources.attempt} {log}
        """

# ---------- preprocess__bowtie2__map ----------
rule preprocess__bowtie2__map:
    """Map one library to reference genome using bowtie2"""
    input:
        forward_=get_input_forward_for_host_mapping,
        reverse_=get_input_reverse_for_host_mapping,
        mock=PRE_INDEX / "{genome}",
        reference=HOSTS / "{genome}.fa.gz",
        faidx=HOSTS / "{genome}.fa.gz.fai",
    output:
        cram=temp(PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram"),
        counts=PRE_QUANT / "{genome}" / "{sample_id}.{library_id}.chr_alignment_counts.tsv"
    log:
        PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.log",
    benchmark:
        PRE_BOWTIE2 / "benchmark/{genome}" / "{sample_id}.{library_id}.tsv"
    params:
        samtools_mem=params["preprocess"]["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    container:
        docker["bowtie2"]
    threads: esc("cpus", "preprocess__bowtie2__map")
    resources:
        runtime=esc("runtime", "preprocess__bowtie2__map"),
        mem_mb=esc("mem_mb", "preprocess__bowtie2__map"),
        cpus_per_task=esc("cpus", "preprocess__bowtie2__map"),
        partition=esc("partition", "preprocess__bowtie2__map"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'preprocess__bowtie2__map')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__bowtie2__map"))
    shell:
        """
        find $(dirname {output.cram}) -name "$(basename {output.cram}).tmp.*.bam" -delete 2> {log}.{resources.attempt} 1>&2

        ( bowtie2 -x {input.mock} -1 {input.forward_} -2 {input.reverse_} \
            --threads {threads} --rg-id '{params.rg_id}' --rg '{params.rg_extra}' \
        | samtools sort -l 9 -M -m {params.samtools_mem} -o {output.cram} \
            --reference {input.reference} --threads {threads} ) \
        2>> {log}.{resources.attempt} 1>&2

        samtools idxstats {output.cram} | awk '{{print $1, $3}}' > {output.counts}

        mv {log}.{resources.attempt} {log}
        """

# ---------- preprocess__bowtie2__extract_nonhost_run ----------
rule preprocess__bowtie2__extract_nonhost_run:
    """Extract non-host reads and convert to FASTQ"""
    input:
        cram=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram",
        reference=HOSTS / "{genome}.fa.gz",
        fai=HOSTS / "{genome}.fa.gz.fai",
    output:
        forward_=temp(PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}_1.fq.gz"),
        reverse_=temp(PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}_2.fq.gz"),
    log:
        PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}.log",
    benchmark:
        PRE_BOWTIE2 / "benchmark/non{genome}" / "{sample_id}.{library_id}.tsv"
    container:
        docker["bowtie2"]
    params:
        samtools_mem=params["preprocess"]["bowtie2"]["samtools"]["mem_per_thread"],
    threads: esc("cpus", "preprocess__bowtie2__extract_nonhost_run")
    resources:
        runtime=esc("runtime", "preprocess__bowtie2__extract_nonhost_run"),
        mem_mb=esc("mem_mb", "preprocess__bowtie2__extract_nonhost_run"),
        cpus_per_task=esc("cpus", "preprocess__bowtie2__extract_nonhost_run"),
        partition=esc("partition", "preprocess__bowtie2__extract_nonhost_run"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'preprocess__bowtie2__extract_nonhost_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__bowtie2__extract_nonhost_run"))
    shell:
        """
        ( samtools view --reference {input.reference} --threads {threads} -u -o /dev/stdout -f 12 {input.cram} \
        | samtools collate -O -u -f --reference {input.reference} -@ {threads} - \
        | samtools fastq -1 >(pigz > {output.forward_}) -2 >(pigz > {output.reverse_}) \
            -0 /dev/null -c 9 --threads {threads} ) 2> {log} 1>&2
        """

# ---------- preprocess__store_final_fastq ----------
rule preprocess__store_final_fastq:
    """Copy final processed FASTQ files"""
    input:
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
    output:
        forward_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
    log:
        PRE_BOWTIE2 / "decontaminated_reads" / "log" / "{sample_id}.{library_id}.log",
    benchmark:
        PRE_BOWTIE2 / "decontaminated_reads" / "benchmark" / "{sample_id}.{library_id}.tsv",
    container:
        docker["bowtie2"]
    threads: esc("cpus", "preprocess__store_final_fastq")
    resources:
        runtime=esc("runtime", "preprocess__store_final_fastq"),
        mem_mb=esc("mem_mb", "preprocess__store_final_fastq"),
        cpus_per_task=esc("cpus", "preprocess__store_final_fastq"),
        partition=esc("partition", "preprocess__store_final_fastq"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'preprocess__store_final_fastq')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__store_final_fastq"))
    shell:
        """
        cp {input.forward_} {output.forward_}
        cp {input.reverse_} {output.reverse_}
        """

# ---------- preprocess__bowtie2__extract_nonhost ----------
rule preprocess__bowtie2__extract_nonhost:
    """Run bowtie2_extract_nonhost for all PE libraries"""
    input:
        [PRE_BOWTIE2 / "decontaminated_reads" / f"{sample_id}.{library_id}_{end}.fq.gz"
         for sample_id, library_id in SAMPLE_LIBRARY
         for end in ["1", "2"]]

# ---------- preprocess__bowtie2 ----------
rule preprocess__bowtie2:
    """Run all the preprocessing steps for bowtie2"""
    input:
        rules.preprocess__bowtie2__extract_nonhost.input
