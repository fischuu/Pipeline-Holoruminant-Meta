rule _annotate__phylophlan_run:
    """
    Run phylophlan for the dereplicated genomes
    """
    input:
        contigs=DREP / "dereplicated_genomes",
        maas=GTDBTK / "gtdbtk.summary.tsv",
        database=lambda w: features["databases"]["phylophlan"][w.phylophlan_db],
    output:
        out_folder=directory(PHYLOPHLAN / "{phylophlan_db}"),
    log:
        PHYLOPHLAN / "{phylophlan_db}.log",
    benchmark:
        PHYLOPHLAN / "benchmark/{phylophlan_db}.tsv",
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        nvme=config["resources"]["nvme"]["large"]
    params:
        diversity=params["annotate"]["phylophlan"]["diversity"],
        config_folder=config["phylophlan-config"],
        out_folder=lambda w: PHYLOPHLAN / w.phylophlan_db,
    conda:
        "__environment__.yml"
    container:
        docker["annotate"]
    shell:
        """
            echo Running Phylophlan on $(hostname) 2>> {log} 1>&2

            mkdir -p {params.out_folder}

            phylophlan -i {input.contigs} \
                       -d {input.database} \
                       -e fa \
                       --diversity {params.diversity} \
                       --nproc {threads} \
                       --output_folder {params.out_folder} \
                       --verbose \
                       2>> {log} 1>&2
        """

rule _annotate__phylophlan:
    """Run phylophlan over all databases."""
    input:
        [
            PHYLOPHLAN / phylophlan_db 
            for phylophlan_db in features["databases"]["phylophlan"]
        ],