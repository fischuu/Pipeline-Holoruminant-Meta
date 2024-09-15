rule _preprocess__diamond__assign:
    """
    Run Diamond over all samples at once using the /dev/shm/ trick.

    NOTE: /dev/shm may be not empty after the job is done.
    """
    input:
        forwards=get_final_forward_from_pre,
        reverses=get_final_reverse_from_pre,
        database=lambda w: features["databases"]["diamond"][w.diamond_db],
    output:
        out_gzs=[
            DIAMOND / "{diamond_db}" / f"{sample_id}.{library_id}.out.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    log:
        DIAMOND / "{diamond_db}.log",
    benchmark:
        DIAMOND / "benchmark/{diamond_db}.tsv",
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"], 
        time =  config["resources"]["time"]["longrun"],
        partition = config["resources"]["partition"]["small"],
        nvme = config["resources"]["nvme"]["large"]
    params:
        in_folder=FASTP,
        out_folder=lambda w: DIAMON / w.diamond_db,
        diamond_db_shm=lambda w:  os.path.join(DIAMONDSHM, w.kraken_db),
    conda:
        "__environment__.yml"
    singularity:
        docker["annotate"]
    shell:
        """
        {{
            echo Running Diamond in $(hostname) 2>> {log} 1>&2

            echo Using quick disc space: {params.diamond_db_shm} 2>> {log} 1>&2

            mkdir --parents {params.diamond_db_shm}
            mkdir --parents {params.out_folder}

            rsync \
                --archive \
                --progress \
                --recursive \
                --times \
                --verbose \
                {input.database}/*.dmnd \
                {params.diamond_db_shm} \
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

rule preprocess__diamond:
    """Run diamond over all samples at once using the /dev/shm/ trick."""
    input:
        [
            DIAMOND / diamond_db / f"{sample_id}.{library_id}.report"
            for sample_id, library_id in SAMPLE_LIBRARY
            for diamond_db in features["databases"]["diamond"]
        ],
