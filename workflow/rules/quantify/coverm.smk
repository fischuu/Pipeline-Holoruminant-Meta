rule quantify__coverm__genome_run:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        cram=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=DREP / "dereplicated_genomes.fa.gz",
        fai=DREP / "dereplicated_genomes.fa.gz.fai",
    output:
        tsv=COVERM / "genome" / "{method}" / "{sample_id}.{library_id}.tsv",
    container:
        docker["quantify"]
    threads: esc("cpus", "quantify__coverm__genome_run")
    resources:
        runtime=esc("runtime", "quantify__coverm__genome_run"),
        mem_mb=esc("mem_mb", "quantify__coverm__genome_run"),
        cpu_per_task=esc("cpus", "quantify__coverm__genome_run"),
        slurm_partition=esc("partition", "quantify__coverm__genome_run"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "quantify__coverm__genome_run", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__coverm__genome_run"))
    log:
        COVERM / "genome" / "{method}" / "{sample_id}.{library_id}.log",
    params:
        method="{method}",
        min_covered_fraction=params["quantify"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["quantify"]["coverm"]["genome"]["separator"],
    shell:
        """
        ( samtools view \
            --exclude-flags 4 \
            --reference {input.reference} \
            --fast \
            {input.cram} \
        | coverm genome \
            --bam-files /dev/stdin \
            --methods {params.method} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output.tsv} \
        ) 2> {log}
        """


rule quantify__coverm__genome_aggregate:
    """Run coverm genome and a single method"""
    input:
        get_tsvs_for_dereplicate_coverm_genome,
    output:
        tsv=COVERM / "genome.{method}.tsv",
    log:
        COVERM / "genome.{method}.log",
    container:
        docker["quantify"]
    params:
        input_dir=lambda w: COVERM / "genome" / w.method,
        script_folder=SCRIPT_FOLDER,
    threads: esc("cpus", "quantify__coverm__genome_aggregate")
    resources:
        runtime=esc("runtime", "quantify__coverm__genome_aggregate"),
        mem_mb=esc("mem_mb", "quantify__coverm__genome_aggregate"),
        cpu_per_task=esc("cpus", "quantify__coverm__genome_aggregate"),
        slurm_partition=esc("partition", "quantify__coverm__genome_aggregate"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "quantify__coverm__genome_aggregate", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__coverm__genome_aggregate"))
    shell:
        """
        Rscript --vanilla {params.script_folder}/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule quantify__coverm__genome:
    """Run coverm genome and all methods"""
    input:
        [
            COVERM / f"genome.{method}.tsv"
            for method in params["quantify"]["coverm"]["genome"]["methods"]
        ],


# coverm contig ----
rule quantify__coverm__contig_one:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        cram=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=QUANT_BOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=DREP / "dereplicated_genomes.fa.gz",
        fai=DREP / "dereplicated_genomes.fa.gz.fai",
    output:
        tsv=COVERM / "contig" / "{method}" / "{sample_id}.{library_id}.tsv",
    container:
        docker["quantify"]
    threads: esc("cpus", "quantify__coverm__contig_one")
    resources:
        runtime=esc("runtime", "quantify__coverm__contig_one"),
        mem_mb=esc("mem_mb", "quantify__coverm__contig_one"),
        cpu_per_task=esc("cpus", "quantify__coverm__contig_one"),
        slurm_partition=esc("partition", "quantify__coverm__contig_one"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "quantify__coverm__contig_one", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__coverm__contig_one"))
    log:
        COVERM / "contig" / "{method}" / "{sample_id}.{library_id}.log",
    params:
        method="{method}",
    shell:
        """
        ( samtools view \
            --exclude-flags 4 \
            --reference {input.reference} \
            --fast \
            {input.cram} \
        | coverm contig \
            --bam-files /dev/stdin \
            --methods {params.method} \
            --proper-pairs-only \
        > {output.tsv} \
        ) 2> {log}
        """


rule quantify__coverm__contig_aggregate:
    """Run coverm contig and a single method"""
    input:
        get_tsvs_for_dereplicate_coverm_contig,
    output:
        tsv=COVERM / "contig.{method}.tsv",
    log:
        COVERM / "contig.{method}.log",
    container:
        docker["quantify"]
    params:
        input_dir=lambda w: COVERM / "contig" / w.method,
        script_folder=SCRIPT_FOLDER,
    threads: esc("cpus", "quantify__coverm__contig_aggregate")
    resources:
        runtime=esc("runtime", "quantify__coverm__contig_aggregate"),
        mem_mb=esc("mem_mb", "quantify__coverm__contig_aggregate"),
        cpu_per_task=esc("cpus", "quantify__coverm__contig_aggregate"),
        slurm_partition=esc("partition", "quantify__coverm__contig_aggregate"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "quantify__coverm__contig_aggregate", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__coverm__contig_aggregate"))
    shell:
        """
        Rscript --vanilla {params.script_folder}/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule quantify__coverm__contig:
    """Run coverm contig and all methods"""
    input:
        [
            COVERM / f"contig.{method}.tsv"
            for method in params["quantify"]["coverm"]["contig"]["methods"]
        ],


rule quantify__coverm:
    input:
        rules.quantify__coverm__contig.input,
        rules.quantify__coverm__genome.input,
