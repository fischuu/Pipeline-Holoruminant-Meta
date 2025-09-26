rule mag_annotate__proteinortho:
    """Run Proteinortho over all bakta output FAA files together"""
    input:
        faa=expand(BAKTAMAG / "bakta_{assembly_id}.faa", assembly_id=ASSEMBLIES),
    output:
        tsv=PROTEINORTHO / "proteinortho.proteinortho.tsv",
    log:
        PROTEINORTHO / "proteinortho.log",
    container:
        docker["mag_annotate"]
    threads: esc("cpus", "mag_annotate__proteinortho")
    resources:
        runtime=esc("runtime", "mag_annotate__proteinortho"),
        mem_mb=esc("mem_mb", "mag_annotate__proteinortho"),
        cpus_per_task=esc("cpus", "mag_annotate__proteinortho"),
        slurm_partition=esc("partition", "mag_annotate__proteinortho"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'mag_annotate__proteinortho')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("mag_annotate__proteinortho"))
    params:
        outdir=PROTEINORTHO
    shell:
        """
        
        mkdir -p {params.outdir}
        
        proteinortho {input.faa} \
            -cpus={threads} \
            -project=proteinortho \
            2>> {log} 1>&2
            
        mv proteinortho.proteinortho.tsv {params.outdir} 
        mv proteinortho.blast-graph {params.outdir}
        mv proteinortho.info {params.outdir}
        mv proteinortho.proteinortho-graph {params.outdir}
        """
