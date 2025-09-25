rule annotate__proteinortho:
    """Run Proteinortho over all bakta output FAA files together"""
    input:
        faa=expand(BAKTAMAG / "bakta_{assembly_id}.faa", assembly_id=ASSEMBLIES),
    output:
        tsv=PROTEINORTHO / "proteinortho.proteinortho.tsv",
    log:
        PROTEINORTHO / "proteinortho.log",
    container:
        docker["annotate"]
    threads: esc("cpus", "annotate__proteinortho")
    resources:
        runtime=esc("runtime", "annotate__proteinortho"),
        mem_mb=esc("mem_mb", "annotate__proteinortho"),
        cpus_per_task=esc("cpus", "annotate__proteinortho"),
        slurm_partition=esc("partition", "annotate__proteinortho"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'annotate__proteinortho')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("annotate__proteinortho"))
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
