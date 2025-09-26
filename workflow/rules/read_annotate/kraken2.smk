rule read_annotate__kraken2__assign:
    """
    Run kraken2 over all samples at once using the /dev/shm/ trick.
    NOTE: /dev/shm may be not empty after the job is done.
    """
    input:
        forwards=[
            FASTP / f"{sample_id}.{library_id}_1.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        reverses=[
            FASTP / f"{sample_id}.{library_id}_2.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        database=lambda w: features["databases"]["kraken2"][w.kraken_db],
    output:
        out_gzs=[
            KRAKEN2 / "{kraken_db}" / f"{sample_id}.{library_id}.out.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        reports=[
            KRAKEN2 / "{kraken_db}" / f"{sample_id}.{library_id}.report"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    log:
        KRAKEN2 / "{kraken_db}.log",
    benchmark:
        KRAKEN2 / "benchmark/{kraken_db}.tsv",
    threads: esc("cpus", "read_annotate__kraken2__assign"),
    resources:
        runtime=esc("runtime", "read_annotate__kraken2__assign"),
        mem_mb=esc("mem_mb", "read_annotate__kraken2__assign"),
        cpus_per_task=esc("cpus", "read_annotate__kraken2__assign"),
        partition=esc("partition", "read_annotate__kraken2__assign"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'read_annotate__kraken2__assign')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("read_annotate__kraken2__assign")),
    params:
        retries=len(get_escalation_order("read_annotate__kraken2__assign")),
        nvme=esc_val("nvme", "read_annotate__kraken2__assign"),
        in_folder=FASTP,
        out_folder=lambda w: KRAKEN2 / w.kraken_db,
        kraken2_shm=KRAKEN2SHM,
        kraken2_nvme=KRAKEN2NVME,
        kraken2_db=lambda w: w.kraken_db,
        kraken_db_shm=lambda w: os.path.join(KRAKEN2SHM, w.kraken_db),
        kraken_db_nvme=lambda w: os.path.join(KRAKEN2NVME, w.kraken_db),
    container:
        docker["preprocess"],
    shell: """
        
        # Don't use the strict mode, to handle nvme assignments in a more safer way
        set +u

        echo "Starting Kraken2 rule attempt {resources.attempt}" > {log}.{resources.attempt} 1>&2
        
        mkdir --parents {params.out_folder} 2>> {log}.{resources.attempt} 1>&2

        DB_SRC="{input.database}"
        SHM="{params.kraken2_shm}"
        DB_SHM="{params.kraken_db_shm}"
        DB_NVME_2={params.kraken2_nvme}/{params.kraken2_db}
        DB_NVME={params.kraken_db_nvme}
        
        : "${{DB_SRC:=}}"
        : "${{DB_SHM:=}}"
        : "${{DB_NVME:=}}"

        echo "DB_SRC: $DB_SRC, DB_SHM: $DB_SHM, DB_NVME: $DB_NVME" 2>> {log}.{resources.attempt} 1>&2

        DB_SIZE=$(du -sb $DB_SRC | cut -f1)
        SHM_AVAIL=$(timeout 10s df --output=avail -B1 "$SHM" 2>/dev/null | tail -1 || echo 0)
        #NVME_AVAIL=$(( {params.nvme} * 2**30 ))
        NVME_AVAIL=$(timeout 10s df --output=avail -B1 "$DB_NVME" 2>/dev/null | tail -1 || echo 0)
        echo "DB_SIZE: $DB_SIZE, SHM_AVAIL: $SHM_AVAIL, NVME_AVAIL: $NVME_AVAIL" 2>> {log}.{resources.attempt} 1>&2

        if [ "$DB_SIZE" -lt "$SHM_AVAIL" ]; then
            DB_DST="$DB_SHM"
        else
            if [ "$DB_SIZE" -lt "$NVME_AVAIL" ]; then
                DB_DST="$DB_NVME"
            else
                if [ {resources.attempt} -eq {params.retries} ]; then
                    DB_DST="{input.database}"
                else
                    echo "DB too large for available storage, aborting attempt {resources.attempt}" 2>> {log}.{resources.attempt} 1>&2
                    exit 1
                fi
            fi
        fi

        if [ "$DB_DST" != "{input.database}" ]; then
            mkdir -p $DB_DST
            rsync --archive --progress --recursive --times --verbose {input.database}/*.k2d $DB_DST 2>> {log}.{resources.attempt}
        fi

        for file in {input.forwards}; do
            sample_id=$(basename $file _1.fq.gz)
            forward={params.in_folder}/${{sample_id}}_1.fq.gz
            reverse={params.in_folder}/${{sample_id}}_2.fq.gz
            output={params.out_folder}/${{sample_id}}.out.gz
            report={params.out_folder}/${{sample_id}}.report
            log_file={params.out_folder}/${{sample_id}}.log

            echo "Processing $sample_id" 2>> {log}.{resources.attempt} 1>&2

            kraken2 \
                --db $DB_DST \
                --threads {threads} \
                --gzip-compressed \
                --paired \
                --output >(pigz --processes {threads} > $output) \
                --report $report \
                --memory-mapping \
                $forward $reverse \
            2> $log_file 1>&2
        done

        if [ "$DB_DST" = "$DB_SHM" ]; then
            rm -rfv $DB_DST 2>> {log}.{resources.attempt}
        fi

        mv {log}.{resources.attempt} {log}
    """


rule read_annotate__kraken2:
    """Run kraken2 over all samples at once using the /dev/shm/ trick."""
    input:
        [
            KRAKEN2 / kraken_db / f"{sample_id}.{library_id}.report"
            for sample_id, library_id in SAMPLE_LIBRARY
            for kraken_db in features["databases"]["kraken2"]
        ],
