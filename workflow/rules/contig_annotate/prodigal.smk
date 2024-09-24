rule _contigAnnotate__prodigal:
    """
    Predict genes from assemblies (PRODIGAL).
    """
    input:
        assembly=lambda wildcards: (
            MEGAHIT / f"{wildcards.assembly_id}.fa.gz" if config["assembler"] == "megahit" else
            METASPADES / f"{wildcards.assembly_id}.fa.gz"
        ),
    output:
        fa=CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.fa",
        gtf=CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.gtf",
        gtfplain=CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal_plain.gtf",
    log:
        CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.log"
    container:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
    shell:""" 
         prodigal -i <(gunzip -c {input.assembly}) \
                  -o {output.gtf} \
                  -a {output.fa} \
                  -p meta -f gff    
                  
         grep -v '^#' {output.gtf} > {output.gtfplain}
    """
    
checkpoint _contigAnnotate__cut_prodigal:
    """
    Cut prodigal output into smaller contigs (BASH).
    """
    input:
        fa=CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.fa",
    output:
        directory(CONTIG_PRODIGAL / "{assembly_id}/Chunks/")
    params:
        out=lambda wildcards: CONTIG_PRODIGAL / f"{wildcards.assembly_id}/Chunks/prodigal.chunk",
        folder=config["pipeline_folder"],
        split=params["contig_annotate"]["prodigal"]["split"]
    log:
        CONTIG_PRODIGAL / "logs/{assembly_id}_cut_prodigal.log"
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["single_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["lowmem"],
        time =  config["resources"]["time"]["shortrun"],
    threads: config["resources"]["cpu_per_task"]["single_thread"]
    container:
        docker["annotate"]
    shell:"""
       mkdir -p {output}
       {params.folder}/workflow/scripts/cutProdigal.sh {params.split} {params.out} {input} 2>> {log} 1>&2
    """    
    
rule contig_assemble__prodigal:
    """Run prodigal on all assemblies"""
    input:
        [CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.fa" for assembly_id in ASSEMBLIES],
