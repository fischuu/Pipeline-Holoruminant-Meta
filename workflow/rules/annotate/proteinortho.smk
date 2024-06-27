rule _annotate__proteinortho_splitfaa:
    """Split bakta output to individual files"""
    input:
        BAKTA / "bakta.faa",
    output:
        directory(BAKTA / "faa"),
    log:
        PROTEINORTHO / "proteinortho.log",
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["single_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["lowmem"],
        time =  config["resources"]["time"]["longrun"],
    params:
        script_folder=config["script_folder"],
    shell:
        """
        {params.script_folder}/split_fasta.sh {input} {output}\
              2>> {log} 1>&2
        """


rule _annotate__proteinortho:
    """Run Proteinortho over the bakta output"""
    input:
        BAKTA / "faa",
    output:
        project=PROTEINORTHO / "proteinortho",
    log:
        PROTEINORTHO / "proteinortho.log",
    conda:
        "annotate.yml"
    singularity:
        docker["annotate"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"]//config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        proteinortho {input}/*.fasta \
                     --cpus {threads} \
                     --project {output.project} \
              2>> {log} 1>&2
        """
