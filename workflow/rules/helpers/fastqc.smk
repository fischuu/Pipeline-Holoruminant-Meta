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
    threads: lambda wc: get_resources(wc, starting_profile="small", escalation_order=["medium"])["cpus"]
    resources:
        time = lambda wc: get_resources(wc, starting_profile="small", escalation_order=["medium"])["time"],
        mem_per_cpu = lambda wc: get_resources(wc, starting_profile="small", escalation_order=["medium"])["mem_mb"],
        cpu_per_task = lambda wc: get_resources(wc, starting_profile="small", escalation_order=["medium"])["cpus"],
        partition = lambda wc: get_resources(wc, starting_profile="small", escalation_order=["medium"])["partition"],
    log:
        "{prefix}_fastqc.log"
    retries: 2  # will escalate from medium → large if first attempt fails
    shell:
        "fastqc {input} 2> {log} 1>&2"
