# For now, I took the GTDBTK annotation out, as we get the taxonomic assignment also from the globlal run

rule annotate__dram_mag__annotate:
    """Annotate dereplicate genomes with DRAM"""
    input:
        contigs=MAGSCOT / "{assembly_id}.fa.gz",
        #gtdbtk_summary=GTDBTK / "gtdbtk.summary.tsv",
        dram_db=features["databases"]["dram"],
    output:
        annotation=DRAMMAG / "{assembly_id}" / "annotate"  / "annotations.tsv",
        trnas=DRAMMAG / "{assembly_id}" / "annotate" / "trnas.tsv",
        rrnas=DRAMMAG / "{assembly_id}" / "annotate" / "rrnas.tsv",
    log:
        DRAM / "{assembly_id}" / "annotate_{assembly_id}.log",
    conda:
        "__environment__.yml"
    container:
        docker["dram"]
    params:
        config=config["dram-config"],
        min_contig_size=1500,
        out_dir=lambda wildcards: f"{DRAMMAG}/{wildcards.assembly_id}",
        tmp_dir=lambda wildcards: f"{DRAMMAG}/{wildcards.assembly_id}/annotate",
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["quitehighmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
        nvme = config["resources"]["nvme"]["small"],
        partition = config["resources"]["partition"]["small"]
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

rule annotate__dram_mag__distill:
    """Distill DRAM annotations."""
    input:
        annotation=DRAMMAG / "{assembly_id}" / "annotate"  / "annotations_{assembly_id}.tsv",
        trnas=DRAMMAG / "{assembly_id}" / "annotate" / "trnas_{assembly_id}.tsv",
        rrnas=DRAMMAG / "{assembly_id}" / "annotate" / "rrnas_{assembly_id}.tsv",
        dram_db=features["databases"]["dram"],
    output:
        genome=DRAMMAG / "{assembly_id}" / "genome_stats.tsv",
        metabolism=DRAMMAG / "{assembly_id}" / "metabolism_summary.xlsx",
        product_html=DRAMMAG / "{assembly_id}" / "product.html",
        product_tsv=DRAMMAG / "{assembly_id}" / "product.tsv",
    log:
        DRAMMAG / "{assembly_id}" / "distill.log2",
    conda:
        "__environment__.yml"
    container:
        docker["dram"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
    params:
        config=config["dram-config"],
        outdir=lambda wildcards: f"{DRAMMAG}/{wildcards.assembly_id}",
        outdir_tmp=lambda wildcards: f"{DRAMMAG}/{wildcards.assembly_id}/distill",
    shell:
        """
        DRAM.py distill \
            --config_loc {params.config} \
            --input_file {input.annotations} \
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
        expand(DRAMMAG / "{assembly_id}" / "annotate"  / "annotations_{assembly_id}.tsv", assembly_id=ASSEMBLIES),
