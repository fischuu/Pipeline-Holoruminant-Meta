rule preprocess__diamond__assign:
    """
    Run Diamond
    """
    input:
        forwards=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverses=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        database=lambda w: features["databases"]["diamond"][w.diamond_db],
    output:
        out_R1= DIAMOND / "{diamond_db}" / "{sample_id}.{library_id}_R1.out",
        out_R2= DIAMOND / "{diamond_db}" / "{sample_id}.{library_id}_R2.out",
    log:
        DIAMOND / "{diamond_db}_{sample_id}_{library_id}.log",
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"], 
        time =  config["resources"]["time"]["longrun"],
        partition = config["resources"]["partition"]["small"],
        nvme = config["resources"]["nvme"]["small"]
    params:
        in_folder=FASTP,
        out_folder=lambda w: DIAMOND / w.diamond_db,
        diamond_db_shm=lambda w:  os.path.join(DIAMONDSHM, w.diamond_db),
        diamond_db_path=lambda w: os.path.join(DIAMONDSHM, w.diamond_db, os.path.basename(features["databases"]["diamond"][w.diamond_db])),
    conda:
        "__environment__.yml"
    container:
        docker["annotate"]
    shell:
        """
            echo Running Diamond in $(hostname) 2>> {log} 1>&2

            echo Using quick disc space: {params.diamond_db_shm} 2>> {log} 1>&2

            mkdir --parents {params.diamond_db_shm}
            mkdir --parents {params.out_folder}

            cp {input.database} {params.diamond_db_shm} 2>> {log} 1>&2

            diamond blastx -d {params.diamond_db_path} \
                           -q {input.forwards} \
                           -o {output.out_R1} \
            2> {log} 1>&2

            diamond blastx -d {params.diamond_db_path} \
                           -q {input.reverses} \
                           -o {output.out_R2} \
           2> {log} 1>&2
        """

rule preprocess__diamond__summarise:
    """Run the appropriate R script after all Diamond processing is complete for a given database."""
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
    threads: config["resources"]["cpu_per_task"]["single_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["single_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
        partition = config["resources"]["partition"]["small"]
    params:
        script=lambda w: config["pipeline_folder"] + "workflow/scripts/downstream_scripts/create_featuretable_diamond_" + w.diamond_db + ".R",
        project_folder=WD,
        database=lambda w: w.diamond_db,
    conda:
        "__environment_r__.yml"
    container:
        docker["r_report"]
    shell:
        """
        R -e "project_folder <- '{params.project_folder}';\
              result_folder <- 'results/preprocess/diamond/{params.database}';\
              source('{params.script}')" &> {log}
        """


rule preprocess__diamond:
    """Run diamond."""
    input:
        [
            DIAMOND / diamond_db / f"{sample_id}.{library_id}_R1.out"
            for sample_id, library_id in SAMPLE_LIBRARY
            for diamond_db in features["databases"]["diamond"]
        ],
        ([
            DIAMOND / diamond_db / "summary_a.tsv"
            for diamond_db in features["databases"]["diamond"]
        ]),
        ([
            DIAMOND / diamond_db / "summary_b.tsv"
            for diamond_db in features["databases"]["diamond"]
        ])

