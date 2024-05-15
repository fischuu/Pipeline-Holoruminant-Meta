rule _preprocess__kraken2__assign:
    """
    Run kraken2 over all samples at once using the /dev/shm/ trick.

    NOTE: /dev/shm may be not empty after the job is done.
    """
    input:
        forwards=[
            FASTP / f"{sample_id}.{library_id}_1.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        rerverses=[
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
    threads: 8
    resources:
        mem_mb=params["preprocess"]["kraken2"]["memory_gb"] * 1024,
        runtime=48 * 60,
    params:
        in_folder=FASTP,
        out_folder=lambda w: KRAKEN2 / w.kraken_db,
        kraken_db_shm="/dev/shm/{kraken_db}",
    conda:
        "__environment__.yml"
    shell:
        """
        {{
            echo Running kraken2 in $(hostname) 2>> {log} 1>&2

            mkdir --parents {params.kraken_db_shm}
            mkdir --parents {params.out_folder}

            rsync \
                --archive \
                --progress \
                --recursive \
                --times \
                --verbose \
                {input.database}/*.k2d \
                {params.kraken_db_shm} \
            2>> {log} 1>&2

            for file in {input.forwards} ; do \

                sample_id=$(basename $file _1.fq.gz)
                forward={params.in_folder}/${{sample_id}}_1.fq.gz
                reverse={params.in_folder}/${{sample_id}}_2.fq.gz
                output={params.out_folder}/${{sample_id}}.out.gz
                report={params.out_folder}/${{sample_id}}.report
                log={params.out_folder}/${{sample_id}}.log

                echo $(date) Processing $sample_id 2>> {log} 1>&2

                kraken2 \
                    --db {params.kraken_db_shm} \
                    --threads {threads} \
                    --gzip-compressed \
                    --paired \
                    --output >(pigz --processes {threads} > $output) \
                    --report $report \
                    --memory-mapping \
                    $forward \
                    $reverse \
                2> $log 1>&2

            done
        }} || {{
            echo "Failed job" 2>> {log} 1>&2
        }}

        rm --force --recursive --verbose {params.kraken_db_shm} 2>>{log} 1>&2
        """


rule preprocess__kraken2:
    """Run kraken2 over all samples at once using the /dev/shm/ trick."""
    input:
        [
            KRAKEN2 / kraken_db / f"{sample_id}.{library_id}.report"
            for sample_id, library_id in SAMPLE_LIBRARY
            for kraken_db in features["databases"]["kraken2"]
        ],
