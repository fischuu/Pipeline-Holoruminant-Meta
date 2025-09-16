ESCALATION = ["small","medium","large"]

rule helpers__fastqc:
    input:
        "{prefix}.fq.gz"
    output:
        html="{prefix}_fastqc.html",
        zip="{prefix}_fastqc.zip"
    container:
        docker["helpers"]
    threads: esc("cpus")
    resources:
        runtime=esc("runtime"),
        mem_mb=esc("mem_mb"),
        cpu_per_task=esc("cpus"),
        partition=esc("partition"),
    log:
        "{prefix}_fastqc.log"
    retries: len(ESCALATION)
    shell:
        "fastqc {input} 2> {log} 1>&2"
