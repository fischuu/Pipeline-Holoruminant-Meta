# ---------- preprocess__diamond__assign ----------
rule preprocess__diamond__assign:
    """Run Diamond"""
    input:
        forwards=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverses=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        database=lambda w: features["databases"]["diamond"][w.diamond_db],
    output:
        out_R1=DIAMOND / "{diamond_db}" / "{sample_id}.{library_id}_R1.out",
        out_R2=DIAMOND / "{diamond_db}" / "{sample_id}.{library_id}_R2.out",
    log:
        DIAMOND / "{diamond_db}_{sample_id}_{library_id}.log",
    threads: esc("cpus", "preprocess__diamond__assign")
    resources:
        runtime=esc("runtime", "preprocess__diamond__assign"),
        mem_mb=esc("mem_mb", "preprocess__diamond__assign"),
        cpu_per_task=esc("cpus", "preprocess__diamond__assign"),
        partition=esc("partition", "preprocess__diamond__assign"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__diamond__assign", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__diamond__assign"))
    params:
        in_folder=FASTP,
        out_folder=lambda w: DIAMOND / w.diamond_db,
        diamond_db_shm=lambda w: os.path.join(DIAMONDSHM, w.diamond_db),
        diamond_db_path=lambda w: os.path.join(DIAMONDSHM, w.diamond_db, os.path.basename(features["databases"]["diamond"][w.diamond_db])),
    container:
        docker["annotate"]
    shell:
        """
        echo Running Diamond in $(hostname) 2>> {log}.{resources.attempt} 1>&2
        echo Using quick disc space: {params.diamond_db_shm} 2>> {log}.{resources.attempt} 1>&2

        mkdir --parents {params.diamond_db_shm}
        mkdir --parents {params.out_folder}

        cp {input.database} {params.diamond_db_shm} 2>> {log}.{resources.attempt} 1>&2

        diamond blastx -d {params.diamond_db_path} -q {input.forwards} -o {output.out_R1} 2> {log}.{resources.attempt} 1>&2
        diamond blastx -d {params.diamond_db_path} -q {input.reverses} -o {output.out_R2} 2> {log}.{resources.attempt} 1>&2

        mv {log}.{resources.attempt} {log}
        """

# ---------- preprocess__diamond__summarise ----------
rule preprocess__diamond__summarise:
    """Run R script to summarise Diamond results"""
    input:
        lambda w: [
            DIAMOND / w.diamond_db / f"{sample_id}.{library_id}_R1.out"
            for sample_id, library_id in SAMPLE_LIBRARY
        ] + [
            DIAMOND / w.diamond_db / f"{sample_id}.{library_id}_R2.out"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
    output:
        DIAMOND / "{diamond_db}" / "summary_a.tsv",
        DIAMOND / "{diamond_db}" / "summary_b.tsv"
    log:
        DIAMOND / "{diamond_db}" / "preprocess__diamond__summarise.log"
    threads: esc("cpus", "preprocess__diamond__summarise")
    resources:
        runtime=esc("runtime", "preprocess__diamond__summarise"),
        mem_mb=esc("mem_mb", "preprocess__diamond__summarise"),
        cpu_per_task=esc("cpus", "preprocess__diamond__summarise"),
        partition=esc("partition", "preprocess__diamond__summarise"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__diamond__summarise", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__diamond__summarise"))
    params:
        script=lambda w: config["pipeline_folder"] + "workflow/scripts/downstream_scripts/create_featuretable_diamond_" + w.diamond_db + ".R",
        project_folder=WD,
        database=lambda w: w.diamond_db,
    container:
        docker["r_report"]
    shell:
        """
        R -e "project_folder <- '{params.project_folder}'; result_folder <- 'results/preprocess/diamond/{params.database}'; source('{params.script}')" &> {log}
        """

# ---------- preprocess__diamond ----------
rule preprocess__diamond:
    """Run all Diamond steps"""
    input:
        [
            DIAMOND / diamond_db / f"{sample_id}.{library_id}_R1.out"
            for sample_id, library_id in SAMPLE_LIBRARY
            for diamond_db in features["databases"]["diamond"]
        ] + [
            DIAMOND / diamond_db / f"{sample_id}.{library_id}_R2.out"
            for sample_id, library_id in SAMPLE_LIBRARY
            for diamond_db in features["databases"]["diamond"]
        ] + [
            DIAMOND / diamond_db / "summary_a.tsv"
            for diamond_db in features["databases"]["diamond"]
        ] + [
            DIAMOND / diamond_db / "summary_b.tsv"
            for diamond_db in features["databases"]["diamond"]
        ]
