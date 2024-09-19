rule _preprocess__diamond__assign:
    """
    Run Diamond
    """
    input:
        forwards=get_final_forward_from_pre,
        reverses=get_final_reverse_from_pre,
        database=lambda w: features["databases"]["diamond"][w.diamond_db],
    output:
        out_R1= DIAMOND / "{diamond_db}" / f"{sample_id}.{library_id}_R1.out",
        out_R2= DIAMOND / "{diamond_db}" / f"{sample_id}.{library_id}_R2.out",
    log:
        DIAMOND / "{diamond_db}.log",
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"], 
        time =  config["resources"]["time"]["longrun"],
        partition = config["resources"]["partition"]["small"],
        nvme = config["resources"]["nvme"]["large"]
    params:
        in_folder=FASTP,
        out_folder=lambda w: DIAMOND / w.diamond_db,
        diamond_db_shm=lambda w:  os.path.join(DIAMONDSHM, w.diamond_db),
    conda:
        "__environment__.yml"
    singularity:
        docker["annotate"]
    shell:
        """
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
                {input.database} \
                {params.diamond_db_shm} \
            2>> {log} 1>&2

            diamond blastx -d {params.diamond_db_shm} \
                           -q {input.forwards} \
                           -o {output.out_R1} \
            2> $log 1>&2

            diamond blastx -d {params.diamond_db_shm} \
                           -q {input.reverses} \
                           -o {output.out_R2} \
           2> $log 1>&2
        """

rule preprocess__diamond:
    """Run diamond over all samples at once using the /dev/shm/ trick."""
    input:
        [
            DIAMOND / "{diamond_db}" / f"{sample_id}.{library_id}_R1.out"
            for sample_id, library_id in SAMPLE_LIBRARY
            for diamond_db in features["databases"]["diamond"]
        ],
