rule _preprocess__krona__visualize:
    """
    Create Krona plots for Kraken2 output.
    """
    input:
        KRONA / "{kraken_db}" / f"{sample_id}.{library_id}.report",
    output:
        KRONA / "{kraken_db}" / f"{sample_id}.{library_id}.html",
    log:
        KRONA / "{kraken_db}.log",
    benchmark:
        KRONA / "benchmark/{kraken_db}.tsv",
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
    conda:
        "__environment__.yml"
    container:
        docker["preprocess"]
    shell:
        """
         ImportTaxomoy.pl -q 2 -t 3 {input} -o {output} 2> $log 1>&2
        """

rule preprocess__krona:
    """Run krona for the kraken output"""
    input:
        [
            KRONA / kraken_db / f"{sample_id}.{library_id}.html"
            for sample_id, library_id in SAMPLE_LIBRARY
            for kraken_db in features["databases"]["kraken2"]
        ],
