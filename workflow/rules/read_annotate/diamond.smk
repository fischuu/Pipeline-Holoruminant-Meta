rule read_annotate__diamond__assign:
    """Run Diamond with SHM/NVME staging if possible."""
    input:
        forwards=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverses=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        database=lambda w: features["databases"]["diamond"][w.diamond_db],
    output:
        out_R1=DIAMOND / "{diamond_db}" / "{sample_id}.{library_id}_R1.out",
        out_R2=DIAMOND / "{diamond_db}" / "{sample_id}.{library_id}_R2.out",
    log:
        DIAMOND / "{diamond_db}_{sample_id}_{library_id}.log",
    benchmark:
        DIAMOND / "benchmark/{diamond_db}_{sample_id}_{library_id}.tsv",
    threads: esc("cpus", "read_annotate__diamond__assign"),
    resources:
        runtime=esc("runtime", "read_annotate__diamond__assign"),
        mem_mb=esc("mem_mb", "read_annotate__diamond__assign"),
        cpus_per_task=esc("cpus", "read_annotate__diamond__assign"),
        partition=esc("partition", "read_annotate__diamond__assign"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'read_annotate__diamond__assign')['nvme']}",
        attempt=lambda wildcards, attempt: attempt,
    retries: len(get_escalation_order("read_annotate__diamond__assign")) - 1,
    params:
        retries=len(get_escalation_order("read_annotate__diamond__assign")) - 1,
        in_folder=FASTP,
        out_folder=lambda w: DIAMOND / w.diamond_db,
        run_in_shm=config["diamond_shm"],
        diamond_shm=DIAMONDSHM,
        diamond_nvme=DIAMONDNVME,
        diamond_db_shm=lambda w: os.path.join(DIAMONDSHM, w.diamond_db),
        diamond_path=lambda w: os.path.basename(features["databases"]["diamond"][w.diamond_db]),
        diamond_db_nvme=lambda w: os.path.join(DIAMONDNVME, w.diamond_db),
    container:
        docker["mag_annotate"],
    shell:
        r"""
        set +u
        echo "Starting Diamond rule attempt {resources.attempt}" > {log}.{resources.attempt} 1>&2

        mkdir --parents {params.out_folder} 2>> {log}.{resources.attempt} 1>&2

        DB_SRC="{input.database}"
        SHM="{params.diamond_shm}"
        NVME="{params.diamond_nvme}"
        DB_SHM="{params.diamond_db_shm}"
        DB_NVME="{params.diamond_db_nvme}"

        if [ "{params.run_in_shm}" != "True" ]; then
            echo "In config/config.yaml is set not to use /dev/shm for this rule" 2>> {log}.{resources.attempt} 1>&2
            SHM=""
            DB_SHM=""
        fi
        

        : "${{DB_SRC:=}}"
        : "${{DB_SHM:=}}"
        : "${{DB_NVME:=}}"
        : "${{SHM:=}}"
        : "${{NVME:=}}"

        set -u

        echo "DB_SRC: $DB_SRC, DB_SHM: $DB_SHM, DB_NVME: $DB_NVME" 2>> {log}.{resources.attempt} 1>&2

        DB_SIZE=$(du -sb $DB_SRC | cut -f1)
        SHM_AVAIL=$(timeout 10s df --output=avail -B1 "$SHM" 2>/dev/null | tail -1 || echo 0)
        NVME_AVAIL=$(timeout 10s df --output=avail -B1 "$NVME" 2>/dev/null | tail -1 || echo 0)

        echo "DB_SIZE: $DB_SIZE, SHM_AVAIL: $SHM_AVAIL, NVME_AVAIL: $NVME_AVAIL" 2>> {log}.{resources.attempt} 1>&2

        if [ "$DB_SIZE" -lt "$SHM_AVAIL" ]; then
            DB_DST="$DB_SHM"/{params.diamond_path}
            echo "Set DB_DST to " $DB_DST 2>> {log}.{resources.attempt} 1>&2
        else
            if [ "$DB_SIZE" -lt "$NVME_AVAIL" ]; then
                DB_DST="$DB_NVME"/{params.diamond_path}
                echo "Set DB_DST to " $DB_DST 2>> {log}.{resources.attempt} 1>&2
            else
                if [ {resources.attempt} -eq {params.retries} ]; then
                    DB_DST="{input.database}"   # use source directory
                    echo "Set DB_DST to " $DB_DST 2>> {log}.{resources.attempt} 1>&2
                else
                    echo "DB too large for available storage, aborting attempt {resources.attempt}" 2>> {log}.{resources.attempt} 1>&2
                    exit 1
                fi
            fi
        fi

        if [ "$DB_DST" != "{input.database}" ]; then
            mkdir -p $DB_DST
            cp $DB_SRC $DB_DST 2>> {log}.{resources.attempt}
            echo "Database copied to DB_DST=$DB_DST" 2>> {log}.{resources.attempt} 1>&2
        fi

        echo "Remove file extension from DB_DST=$DB_DST" 2>> {log}.{resources.attempt} 1>&2
        DB_DST_BASE="${DB_DST%.dmnd}"

        echo "Running Diamond using DB_DST_BASE=$DB_DST_BASE" 2>> {log}.{resources.attempt} 1>&2

        diamond blastx -d $DB_DST_BASE -q {input.forwards} -o {output.out_R1} -p {threads} 2>> {log}.{resources.attempt}
        diamond blastx -d $DB_DST_BASE -q {input.reverses} -o {output.out_R2} -p {threads} 2>> {log}.{resources.attempt}

        if [ "$DB_DST" = "$DB_SHM" ]; then
            rm -rfv $DB_DST 2>> {log}.{resources.attempt}
        fi

        mv {log}.{resources.attempt} {log}
        """


rule read_annotate__diamond__summarise:
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
        DIAMOND / "{diamond_db}" / "read_annotate__diamond__summarise.log"
    threads: esc("cpus", "read_annotate__diamond__summarise")
    resources:
        runtime=esc("runtime", "read_annotate__diamond__summarise"),
        mem_mb=esc("mem_mb", "read_annotate__diamond__summarise"),
        cpus_per_task=esc("cpus", "read_annotate__diamond__summarise"),
        partition=esc("partition", "read_annotate__diamond__summarise"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'read_annotate__diamond__summarise')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("read_annotate__diamond__summarise"))
    params:
        script=lambda w: config["pipeline_folder"] + "workflow/scripts/downstream_scripts/create_featuretable_diamond_" + w.diamond_db + ".R",
        project_folder=WD,
        database=lambda w: w.diamond_db,
    container:
        docker["r_report"]
    shell:
        """
        R -e "project_folder <- '{params.project_folder}'; result_folder <- 'results/read_annotate/diamond/{params.database}'; source('{params.script}')" &> {log}
        """

rule read_annotate__diamond:
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
