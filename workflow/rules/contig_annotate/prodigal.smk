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
    singularity:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
    shell:""" 
         prodigal -i {input.assembly} \
                  -o {output.gtf} \
                  -a {output.fa} \
                  -p meta -f gff    
                  
         sed '/^#/d' {output.gtf} > {output.gtfplain}
    """
    
rule contig_assemble__prodigal:
    """Run prodigal on all assemblies"""
    input:
        [CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.fa" for assembly_id in ASSEMBLIES],
