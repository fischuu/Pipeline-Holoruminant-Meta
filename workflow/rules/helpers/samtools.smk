ESCALATION = ["medium","large"]

rule helpers__samtools__index_bam:
    """Index a bam file"""
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    container:
        docker["helpers"]
    threads: esc("cpus")
    resources:
        runtime=esc("runtime"),
        mem_mb=esc("mem_mb"),
        cpu_per_task=esc("cpus"),
        partition=esc("partition"),
    retries: len(ESCALATION)
    log:
        "{prefix}.bam.bai.log",
    benchmark:
        "benchmark/{prefix}.bam.bai.tsv",
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule helpers__samtools__index_cram:
    """Index a cram file"""
    input:
        "{prefix}.cram",
    output:
        "{prefix}.cram.crai",
    container:
        docker["helpers"]
    threads: esc("cpus")
    resources:
        runtime=esc("runtime"),
        mem_mb=esc("mem_mb"),
        cpu_per_task=esc("cpus"),
        partition=esc("partition"),
    retries: len(ESCALATION)
    log:
        "{prefix}.cram.crai.log",
    benchmark:
        "benchmark/{prefix}.cram.crai.tsv",
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule helpers__samtools__faidx_fa:
    """Index a fasta file"""
    input:
        "{prefix}.fa",
    output:
        "{prefix}.fa.fai",
    container:
        docker["helpers"]
    threads: esc("cpus")
    resources:
        runtime=esc("runtime"),
        mem_mb=esc("mem_mb"),
        cpu_per_task=esc("cpus"),
        partition=esc("partition"),
    retries: len(ESCALATION)
    log:
        "{prefix}.fa.fai.log",
    benchmark:
        "benchmark/{prefix}.fa.fai.tsv",
    shell:
        "samtools faidx {input} 2> {log} 1>&2"


rule helpers__samtools__faidx_fagz:
    """Index a gzipped fasta file"""
    input:
        "{prefix}.fa.gz",
    output:
        fai="{prefix}.fa.gz.fai",
        gzi="{prefix}.fa.gz.gzi",
    container:
        docker["helpers"]
    threads: esc("cpus")
    resources:
        runtime=esc("runtime"),
        mem_mb=esc("mem_mb"),
        cpu_per_task=esc("cpus"),
        partition=esc("partition"),
    retries: len(ESCALATION)
    log:
        "{prefix}.fa.gz.log",
    benchmark:
        "benchmark/{prefix}.fa.gz.tsv",
    shell:
        "samtools faidx {input} 2> {log} 1>&2"


rule helpers__samtools__idxstats_cram:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    benchmark:
        "benchmark/{prefix}.idxstats.tsv",
    container:
        docker["helpers"]
    threads: esc("cpus")
    resources:
        runtime=esc("runtime"),
        mem_mb=esc("mem_mb"),
        cpu_per_task=esc("cpus"),
        partition=esc("partition"),
    retries: len(ESCALATION)
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"

ESCALATION = ["small","medium","large"]

rule helpers__samtools__flagstats_cram:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    container:
        docker["helpers"]
    threads: esc("cpus")
    resources:
        runtime=esc("runtime"),
        mem_mb=esc("mem_mb"),
        cpu_per_task=esc("cpus"),
        partition=esc("partition"),
    retries: len(ESCALATION)
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"
