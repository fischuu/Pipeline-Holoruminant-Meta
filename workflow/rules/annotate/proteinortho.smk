rule annotate__proteinortho:
    """Run Proteinortho over all bakta output FAA files together"""
    input:
        faa=expand(BAKTAMAG / "bakta_{assembly_id}.faa", assembly_id=ASSEMBLIES),
    output:
        project=PROTEINORTHO / "proteinortho.proteinortho.tsv",
    log:
        PROTEINORTHO / "proteinortho.log",
    container:
        docker["annotate"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"]//config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"]
    shell:
        """
        proteinortho {input.faa} \
            -cpus={threads} \
            -project={output.project} \
            2>> {log} 1>&2
        """
