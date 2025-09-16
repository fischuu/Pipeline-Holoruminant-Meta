rule helpers__fastqc:
    input:
        "{prefix}.fq.gz"
    output:
        html="{prefix}_fastqc.html",
        zip="{prefix}_fastqc.zip"
    conda:
        "__environment__.yml"
    container:
        docker["helpers"]
    threads: lambda wc, attempt: get_resources(wc, attempt, escalation_order=["small", "medium"])["cpus"]
    resources:
        time = lambda wc, attempt: get_resources(wc, attempt, escalation_order=["small","medium"])["time"],
        mem_per_cpu = lambda wc, attempt: get_resources(wc, attempt, escalation_order=["small","medium"])["mem_mb"],
        cpu_per_task = lambda wc, attempt: get_resources(wc, attempt, escalation_order=["small","medium"])["cpus"],
        partition = lambda wc, attempt: get_resources(wc, attempt, escalation_order=["small","medium"])["partition"],
    log:
        "{prefix}_fastqc.log"
    retries: 2  # will escalate two positions through the escalation order (meaning, set this number to the length of above escalation order...)
    shell:
        "fastqc {input} 2> {log} 1>&2"
