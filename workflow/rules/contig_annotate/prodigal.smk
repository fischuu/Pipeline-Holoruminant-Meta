rule contig_annotate__prodigal_run:
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
    threads: esc("cpus", "contig_annotate__prodigal_run")
    resources:
        runtime=esc("runtime", "contig_annotate__prodigal_run"),
        mem_mb=esc("mem_mb", "contig_annotate__prodigal_run"),
        cpus_per_task=esc("cpus", "contig_annotate__prodigal_run"),
        slurm_partition=esc("partition", "contig_annotate__prodigal_run"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'contig_annotate__prodigal_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("contig_annotate__prodigal_run"))
    params:
        tmp_file=lambda wildcards: f"{CONTIG_PRODIGAL}/{wildcards.assembly_id}.fa",
    shell:""" 
         prodigal -i <(gunzip -c {input.assembly}) \
                  -o {output.gtf} \
                  -a {output.fa} \
                  -p meta -f gff
                  
         grep -v '^#' {output.gtf} > {output.gtfplain}
    """
    
checkpoint contig_annotate__cut_prodigal:
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
    threads: esc("cpus", "contig_annotate__cut_prodigal")
    resources:
        runtime=esc("runtime", "contig_annotate__cut_prodigal"),
        mem_mb=esc("mem_mb", "contig_annotate__cut_prodigal"),
        cpus_per_task=esc("cpus", "contig_annotate__cut_prodigal"),
        slurm_partition=esc("partition", "contig_annotate__cut_prodigal"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'contig_annotate__cut_prodigal')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("contig_annotate__cut_prodigal"))
    container:
        docker["mag_annotate"]
    shell:"""
       mkdir -p {output}
       {params.folder}/workflow/scripts/cutProdigal.sh {params.split} {params.out} {input} 2>> {log} 1>&2
    """    
    
rule contig_annotate__prodigal:
    """Run prodigal on all assemblies"""
    input:
       # [CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.fa" for assembly_id in ASSEMBLIES],
        expand(CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.fa", assembly_id=ASSEMBLIES),
