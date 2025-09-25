rule report__step__reads:
    """Collect all reports for the reads step"""
    input:
        rules.reads__fastqc.input,
    output:
        html=REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    container:
        docker["report"]
    params:
        dir=REPORT_STEP,
    threads: esc("cpus", "report__step__reads")
    resources:
        runtime=esc("runtime", "report__step__reads"),
        mem_mb=esc("mem_mb", "report__step__reads"),
        cpus_per_task=esc("cpus", "report__step__reads"),
        slurm_partition=esc("partition", "report__step__reads"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "report__step__reads", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__step__reads"))
    shell:
        """
        multiqc \
            --filename reads \
            --title reads \
            --force \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__preprocess:
    """Collect all reports for the preprocessing step"""
    input:
        rules.preprocess__fastp.input.json,
        rules.preprocess__fastqc.input,
        rules.preprocess__samtools.input,
        rules.preprocess__kraken2.input,
    output:
        html=REPORT_STEP / "preprocess.html",
    log:
        REPORT_STEP / "preprocess.log",
    container:
        docker["report"]
    params:
        dir=REPORT_STEP,
    threads: esc("cpus", "report__step__preprocess")
    resources:
        runtime=esc("runtime", "report__step__preprocess"),
        mem_mb=esc("mem_mb", "report__step__preprocess"),
        cpus_per_task=esc("cpus", "report__step__preprocess"),
        slurm_partition=esc("partition", "report__step__preprocess"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "report__step__preprocess", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__step__preprocess"))
    shell:
        """
        multiqc \
            --title preprocess \
            --force \
            --filename preprocess \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input} \
        2> {log}.{resources.attempt} 1>&2

        mv {log}.{resources.attempt} {log}
        """


rule report__step__assemble:
    """Collect all reports from the assemble step"""
    input:
        QUAST,
    output:
        REPORT_STEP / "assemble.html",
    log:
        REPORT_STEP / "assemble.log",
    container:
        docker["report"]
    params:
        dir=REPORT_STEP,
    threads: esc("cpus", "report__step__assemble")
    resources:
        runtime=esc("runtime", "report__step__assemble"),
        mem_mb=esc("mem_mb", "report__step__assemble"),
        cpus_per_task=esc("cpus", "report__step__assemble"),
        slurm_partition=esc("partition", "report__step__assemble"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "report__step__assemble", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__step__assemble"))
    shell:
        """
        multiqc \
            --title assemble \
            --force \
            --filename assemble \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__quantify:
    """Collect all reports from the quantify step"""
    input:
        rules.quantify__samtools.input,
    output:
        REPORT_STEP / "quantify.html",
    log:
        REPORT_STEP / "quantify.log",
    container:
        docker["report"]
    params:
        dir=REPORT_STEP,
    threads: esc("cpus", "report__step__quantify")
    resources:
        runtime=esc("runtime", "report__step__quantify"),
        mem_mb=esc("mem_mb", "report__step__quantify"),
        cpus_per_task=esc("cpus", "report__step__quantify"),
        slurm_partition=esc("partition", "report__step__quantify"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "report__step__quantify", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__step__quantify"))
    shell:
        """
        multiqc \
            --title quantify \
            --force \
            --filename quantify \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step:
    """Report for all steps"""
    input:
        REPORT_STEP / "reads.html",
        REPORT_STEP / "preprocess.html",
        REPORT_STEP / "assemble.html",
        REPORT_STEP / "quantify.html",
