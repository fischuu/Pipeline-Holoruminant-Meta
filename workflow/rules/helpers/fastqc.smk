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
        cpus_per_task=esc("cpus", "helpers__fastqc"),
        slurm_partition=esc("partition", "helpers__fastqc"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'helpers__fastqc')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("helpers__fastqc"))
    log:
        "{prefix}_fastqc.log"
    shell:
        "fastqc {input} 2> {log} 1>&2"
