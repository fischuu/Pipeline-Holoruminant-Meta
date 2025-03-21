rule preprocess__krona__visualize:
    """
    Create Krona plots for Kraken2 output.
    """
    input:
        KRAKEN2 / "{kraken_db}" / "{sample_id}.{library_id}.report",
    output:
        KRONA / "{kraken_db}" / "{sample_id}.{library_id}.html",
    log:
        KRONA / "{kraken_db}_{sample_id}.{library_id}.log",
    benchmark:
        KRONA / "benchmark/{kraken_db}_{sample_id}.{library_id}bsk.tsv",
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
    conda:
        "__environment__.yml"
    container:
        docker["preprocess"]
    params:
        krona_db = features["databases"]["krona"],
    shell:
        """
         ktImportTaxonomy -q 2 -t 3 {input} -o {output} --tax {params.krona_db} 2> {log} 1>&2
        """

rule preprocess__krona:
    """Run krona for the kraken output"""
    input:
        [
            KRONA / kraken_db / f"{sample_id}.{library_id}.html"
            for sample_id, library_id in SAMPLE_LIBRARY
            for kraken_db in features["databases"]["kraken2"]
        ],
