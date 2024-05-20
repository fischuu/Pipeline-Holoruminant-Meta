rule _helpers__fastqc:
    """Run fastqc on a single file"""
    input:
        "{prefix}.fq.gz",
    output:
        html="{prefix}_fastqc.html",
        zip="{prefix}_fastqc.zip",
    conda:
        "__environment__.yml"
    singularity:
        docker["helpers"]
    resources:
        time =  config["resources"]["time"]["shortrun"]
    log:
        "{prefix}_fastqc.log",
    shell:
        "fastqc {input} 2> {log} 1>&2"
