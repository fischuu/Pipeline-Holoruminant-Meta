rule preprocess__metaphlan__run:
    """Run metaphlan over one sample
    """
    input:
        forward_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        data=features["databases"]["metaphlan4"],
    output:
        bt2_out=METAPHLAN / "bowtie2_out" / "{sample_id}.{library_id}.bz2",
        mp_out=METAPHLAN / "profiled" / "{sample_id}.{library_id}.txt",
    log:
        METAPHLAN / "log" / "{sample_id}.{library_id}.log",
    benchmark:
        METAPHLAN / "benchmark" / "{sample_id}.{library_id}.tsv",
    container:
        docker["preprocess"]
    threads: esc("cpus", "preprocess__metaphlan__run")
    resources:
        runtime=esc("runtime", "preprocess__metaphlan__run"),
        mem_mb=esc("mem_mb", "preprocess__metaphlan__run"),
        cpus_per_task=esc("cpus", "preprocess__metaphlan__run"),
        slurm_partition=esc("partition", "preprocess__metaphlan__run"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__metaphlan__run", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__metaphlan__run"))
    params:
        tmp = config["tmp_storage"]
    shell:
        """
        TMPDIR={params.tmp}

        metaphlan {input.forward_},{input.reverse_} \
                  --bowtie2out {output.bt2_out} \
                  --nproc {threads} \
                  --input_type fastq \
                  -o {output.mp_out} \
                  --bowtie2db {input.data} \
        2>> {log} 1>&2
        """


rule preprocess__metaphlan__condense:
    """Aggregate all the metaphlan results into a single table"""
    input:
        profiled_data=[
            METAPHLAN / "profiled" / f"{sample_id}.{library_id}.txt"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
    output:
        METAPHLAN / "metaphlan_profiled.tsv",
    log:
        METAPHLAN / "metaphlan.log",
    benchmark:
        METAPHLAN / "benchmark/metaphlan.tsv",
    container:
        docker["preprocess"]
    params:
        input_dir=METAPHLAN,
    threads: esc("cpus", "preprocess__metaphlan__condense")
    resources:
        runtime=esc("runtime", "preprocess__metaphlan__condense"),
        mem_mb=esc("mem_mb", "preprocess__metaphlan__condense"),
        cpus_per_task=esc("cpus", "preprocess__metaphlan__condense"),
        slurm_partition=esc("partition", "preprocess__metaphlan__condense"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__metaphlan__condense", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__metaphlan__condense"))
    shell:
        """
         merge_metaphlan_tables.py {input.profiled_data} > {output}
        #cat {input.profiled_data} > {output}
        #2> {log} 1>&2
        """

rule preprocess__metaphlan:
    input:
        rules.preprocess__metaphlan__condense.output
