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
    threads: get_resources(escalation_order=["small", "medium", "large"])["cpus"]
    resources:
        runtime = get_resources(escalation_order=["small","medium","large"])["runtime"],
        mem_mb = get_resources(escalation_order=["small","medium","large"])["mem_mb"],
        cpu_per_task = get_resources(escalation_order=["small","medium","large"])["cpus"],
        partition = get_resources(escalation_order=["small","medium","large"])["partition"],
    log:
        "{prefix}_fastqc.log"
    retries: 3  # will escalate three positions through the escalation order (meaning, set this number to the length of above escalation order...)
    shell:
        "fastqc {input} 2> {log} 1>&2"
