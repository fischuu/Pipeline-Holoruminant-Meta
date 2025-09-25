# For now, I took the GTDBTK annotation out, as we get the taxonomic assignment also from the globlal run

rule annotate__dram_mag__annotate:
    """Annotate dereplicate genomes with DRAM"""
    input:
        contigs=MAGSCOT / "{assembly_id}.fa.gz",
        #gtdbtk_summary=GTDBTK / "gtdbtk.summary.tsv",
        dram_db=features["databases"]["dram"],
    output:
        annotation=DRAMMAG / "{assembly_id}" / "annotate"  / "{assembly_id}_annotations.tsv",
        trnas=DRAMMAG / "{assembly_id}" / "annotate" / "{assembly_id}_trnas.tsv",
        rrnas=DRAMMAG / "{assembly_id}" / "annotate" / "{assembly_id}_rrnas.tsv",
    log:
        DRAM / "{assembly_id}" / "annotate_{assembly_id}.log",
    container:
        docker["dram"]
    params:
        config=config["dram-config"],
        min_contig_size=1500,
        out_dir=lambda wildcards: f"{DRAMMAG}/{wildcards.assembly_id}",
        tmp_dir=lambda wildcards: f"{DRAMMAG}/{wildcards.assembly_id}/annotate",
    threads: esc("cpus", "annotate__dram_mag__annotate")
    resources:
        runtime=esc("runtime", "annotate__dram_mag__annotate"),
        mem_mb=esc("mem_mb", "annotate__dram_mag__annotate"),
        cpus_per_task=esc("cpus", "annotate__dram_mag__annotate"),
        slurm_partition=esc("partition", "annotate__dram_mag__annotate"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "annotate__dram_mag__annotate", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("annotate__dram_mag__annotate"))
    shell:
        """
        rm -rf {params.tmp_dir}
        
        echo "Hostname: $(hostname)" 2>> {log} 1>&2
        echo "Temporary directory: $TMPDIR" 2>> {log} 1>&2
        df -h 2>> {log} 1>&2
        
        DRAM.py annotate \
                --config_loc {params.config} \
                --input_fasta {input.contigs} \
                --output_dir {params.tmp_dir} \
                --threads {threads} \
        2>> {log} 1>&2
    """
    
rule annotate__fix_dram_mag_annotations_scaffold:
    input:
        DRAMMAG / "{assembly_id}" / "annotate"  / "{assembly_id}_annotations.tsv",
    output:
        DRAMMAG / "{assembly_id}" / "annotate"  / "{assembly_id}_annotations.fixed.tsv",
    log:
        DRAMMAG / "{assembly_id}_fix_annotation.log",
    container:
        docker["assemble"]
    params:
        script_folder=SCRIPT_FOLDER,
    shell:
        """
        python {params.script_folder}/fix_annotations.py {input} {output} 2>> {log} 1>&2
        """
        
    

rule annotate__dram_mag__distill:
    """Distill DRAM annotations."""
    input:
        annotation=DRAMMAG / "{assembly_id}" / "annotate"  / "{assembly_id}_annotations.fixed.tsv",
        trnas=DRAMMAG / "{assembly_id}" / "annotate" / "{assembly_id}_trnas.tsv",
        rrnas=DRAMMAG / "{assembly_id}" / "annotate" / "{assembly_id}_rrnas.tsv",
        dram_db=features["databases"]["dram"],
    output:
        genome=DRAMMAG / "{assembly_id}" / "genome_stats.tsv",
        metabolism=DRAMMAG / "{assembly_id}" / "metabolism_summary.xlsx",
        product_html=DRAMMAG / "{assembly_id}" / "product.html",
        product_tsv=DRAMMAG / "{assembly_id}" / "product.tsv",
    log:
        DRAMMAG / "{assembly_id}" / "distill.log2",
    container:
        docker["dram"]
    threads: esc("cpus", "annotate__fix_dram_mag_annotations_scaffold")
    resources:
        runtime=esc("runtime", "annotate__fix_dram_mag_annotations_scaffold"),
        mem_mb=esc("mem_mb", "annotate__fix_dram_mag_annotations_scaffold"),
        cpus_per_task=esc("cpus", "annotate__fix_dram_mag_annotations_scaffold"),
        slurm_partition=esc("partition", "annotate__fix_dram_mag_annotations_scaffold"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "annotate__fix_dram_mag_annotations_scaffold", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("annotate__fix_dram_mag_annotations_scaffold"))
    params:
        config=config["dram-config"],
        outdir=lambda wildcards: f"{DRAMMAG}/{wildcards.assembly_id}",
        outdir_tmp=lambda wildcards: f"{DRAMMAG}/{wildcards.assembly_id}/distill",
    shell:
        """
        DRAM.py distill \
            --config_loc {params.config} \
            --input_file {input.annotation} \
            --rrna_path {input.rrnas} \
            --trna_path {input.trnas} \
            --output_dir {params.outdir_tmp} \
        2> {log} 1>&2

        mv {params.outdir_tmp}/* {params.outdir}/ 2>> {log} 1>&2
        rmdir {params.outdir_tmp} 2>> {log} 1>&2
        """

rule annotate__dram_mags:
    """Run Bakta over the dereplicated mags"""
     input:
        #expand(DRAMMAG / "{assembly_id}" / "annotate"  / "{assembly_id}_annotations.tsv", assembly_id=ASSEMBLIES),
        expand(DRAMMAG / "{assembly_id}" / "genome_stats.tsv", assembly_id=ASSEMBLIES),
