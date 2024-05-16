rule _helpers__samtools__index_bam:
    """Index a bam file"""
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    conda:
        "__environment__.yml"
    singularity:
        docker["helpers"]
    log:
        "{prefix}.bam.bai.log",
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule _helpers__samtools__index_cram:
    """Index a cram file"""
    input:
        "{prefix}.cram",
    output:
        "{prefix}.cram.crai",
    conda:
        "__environment__.yml"
    singularity:
        docker["helpers"]
    log:
        "{prefix}.cram.crai.log",
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule _helpers__samtools__faidx_fa:
    """Index a fasta file"""
    input:
        "{prefix}.fa",
    output:
        "{prefix}.fa.fai",
    conda:
        "__environment__.yml"
    singularity:
        docker["helpers"]
    log:
        "{prefix}.fa.fai.log",
    shell:
        "samtools faidx {input} 2> {log} 1>&2"


rule _helpers__samtools__faidx_fagz:
    """Index a gzipped fasta file"""
    input:
        "{prefix}.fa.gz",
    output:
        fai="{prefix}.fa.gz.fai",
        gzi="{prefix}.fa.gz.gzi",
    conda:
        "__environment__.yml"
    singularity:
        docker["helpers"]
    log:
        "{prefix}.fa.gz.log",
    shell:
        "samtools faidx {input} 2> {log} 1>&2"


rule _helpers__samtools__idxstats_cram:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["helpers"]
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"


rule _helpers__samtools__flagstats_cram:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["helpers"]
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"
