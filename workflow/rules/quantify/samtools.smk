rule _quantify__samtools__stats_cram:
    """Get stats from CRAM files using samtools stats."""
    input:
        cram=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=DREP / "dereplicated_genomes.fa.gz",
        fai=DREP / "dereplicated_genomes.fa.gz.fai",
    output:
        txt=QUANT_BOWTIE2 / "{sample_id}.{library_id}.stats.txt",
    log:
        QUANT_BOWTIE2 / "{sample_id}.{library_id}.stats.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["quantify"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        "samtools stats --reference {input.reference} {input.cram} > {output.txt} 2> {log}"


rule quantify__samtools:
    """Get stats from CRAM files using samtools stats."""
    input:
        [
            QUANT_BOWTIE2 / f"{sample_id}.{library_id}.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for extension in ["stats.txt", "flagstats.txt"]
        ],
