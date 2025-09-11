rule contig_annotate__sylph_profile:
    """Run Sylph profiler"""
    input:
        contigs=CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.fa",
        db=features["databases"]["sylph"],
    output:
        CONTIG_SYLPH / "{assembly_id}" / "profiling.tsv",
    log:
        CONTIG_SYLPH / "{assembly_id}" / "quality_report.log",
    conda:
        "sylph.yml"
    container:
        docker["sylph"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"]//config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
        nvme = config["resources"]["nvme"]["small"],
    shell:
        """
        sylph profile {input.db} -1 {input.contigs} -t {threads} > {output} 2>> {log} 1>&2
        """


rule contig_annotate__sylph:
    """Run Sylph"""
    input:
        [CONTIG_SYLPH / "{assembly_id}" / "profiling.tsv" for assembly_id in ASSEMBLIES],
