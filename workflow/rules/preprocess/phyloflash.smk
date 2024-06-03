rule _preprocess__PhyloFlash__run:
    """Run PhyloFlash over one sample
    """
    input:
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        phyloflash_dbs=PHYLOFLASH_DBS,
    output:
        phyloflash_out=PHYLOFLASH / "{sample_id}.{library_id}_genefamilies.tsv",
    log:
        PHYLOFLASH / "log" / "{sample_id}.{library_id}.log",
    benchmark:
        PHYLOFLASH / "benchmark" / "{sample_id}.{library_id}.tsv",
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
    threads: config["resources"]["mem_per_cpu"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["mem_per_cpu"]["multi_thread"]
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["mem_per_cpu"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"]
    shell:
        """
        phyloFlash.pl -dbhome {input.phyloflash_dbs} \
                      -lib Phyloflash -CPUs {threads} \
                      -read1 {input.forward_} \
                      -read2 {input.reverse_} \
                      -almosteverything
        """


rule _preprocess__phyloflash__condense:
    """Aggregate all the PhyloFlash results into a single table"""
    input:
        genefamily_data=[
            PHYLOFLASH / f"{sample_id}.{library_id}_genefamilies.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
    output:
        PHYLOFLASH / "phyloflash_genefamilies.tsv",
    log:
        PHYLOFLASH / "phyloflash.log",
    benchmark:
        PHYLOFLASH / "benchmark/phyloflash.tsv",
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
    params:
        input_dir=PHYLOFLASH,
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"]
    shell:
        """
        cat {input.genefamily_data} > {output}
        2> {log} 1>&2
        """

rule preprocess__phyloflash:
    input:
        rules._preprocess__phyloflash__condense.output
