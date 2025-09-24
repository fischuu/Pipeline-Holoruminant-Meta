rule helpers__fastqc:
    input:
        "{prefix}.fq.gz"
    output:
        html="{prefix}_fastqc.html",
        zip="{prefix}_fastqc.zip"
    container:
        docker["helpers"]
    threads: esc("cpus", "helpers__fastqc")
    resources:
        runtime=esc("runtime", "helpers__fastqc"),
        mem_mb=esc("mem_mb", "helpers__fastqc"),
        cpu_per_task=esc("cpus", "helpers__fastqc"),
        slurm_partition=esc("partition", "helpers__fastqc"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "helpers__fastqc", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("helpers__fastqc"))
    log:
        "{prefix}_fastqc.log"
    shell:
        "fastqc {input} 2> {log} 1>&2"
