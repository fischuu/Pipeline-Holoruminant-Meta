rule preprocess__phyloflash__run:
    """Run PhyloFlash over one sample
    """
    input:
        forward_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        phyloflash_dbs=PHYLOFLASH_DBS,
    output:
        phyloflash_out=PHYLOFLASH / "{sample_id}_{library_id}.phyloFlash.html",
        phyloflash_log=PHYLOFLASH / "{sample_id}_{library_id}.phyloFlash.log",
    log:
        PHYLOFLASH / "log" / "{sample_id}.{library_id}.log",
    benchmark:
        PHYLOFLASH / "benchmark" / "{sample_id}.{library_id}.tsv",
    container:
        docker["phyloflash"]
    threads: esc("cpus", "preprocess__phyloflash__run")
    resources:
        runtime=esc("runtime", "preprocess__phyloflash__run"),
        mem_mb=esc("mem_mb", "preprocess__phyloflash__run"),
        cpus_per_task=esc("cpus", "preprocess__phyloflash__run"),
        slurm_partition=esc("partition", "preprocess__phyloflash__run"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__phyloflash__run", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__phyloflash__run"))
    params:
        lib="{sample_id}_{library_id}",
        outdir=PHYLOFLASH
    shell:
        """
        echo "Using temporary directory: $TMPDIR" 2> {log} 1>&2
        df -h $TMPDIR 2>> {log} 1>&2
        ls -lah $TMPDIR 2>> {log} 1>&2

        
        phyloFlash.pl -dbhome {input.phyloflash_dbs} \
                      -lib {params.lib} -CPUs {threads} \
                      -read1 {input.forward_} \
                      -read2 {input.reverse_} \
                      -almosteverything 2>> {log} 1>&2

        mkdir -p {params.outdir}
        mv {params.lib}.* {params.outdir}
        """


rule preprocess__phyloflash__condense:
    """Aggregate all the PhyloFlash results into a single table"""
    input:
        genefamily_data=[
            PHYLOFLASH / f"{sample_id}_{library_id}.phyloFlash.log"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
    output:
        PHYLOFLASH / "phyloflash_all.log",
    log:
        PHYLOFLASH / "phyloflash.log",
    benchmark:
        PHYLOFLASH / "benchmark/phyloflash.tsv",
    container:
        docker["phyloflash"]
    params:
        input_dir=PHYLOFLASH,
    threads: esc("cpus", "preprocess__phyloflash__condense")
    resources:
        runtime=esc("runtime", "preprocess__phyloflash__condense"),
        mem_mb=esc("mem_mb", "preprocess__phyloflash__condense"),
        cpus_per_task=esc("cpus", "preprocess__phyloflash__condense"),
        slurm_partition=esc("partition", "preprocess__phyloflash__condense"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__phyloflash__condense", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__phyloflash__condense"))
    shell:
        """
        cat {input.genefamily_data} > {output}
        2> {log} 1>&2
        """

rule preprocess__phyloflash:
    input:
        rules.preprocess__phyloflash__condense.output
