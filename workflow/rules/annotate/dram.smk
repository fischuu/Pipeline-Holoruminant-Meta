rule annotate__dram__annotate:
    """Annotate dereplicate genomes with DRAM"""
    input:
        dereplicated_genomes=DREP / "dereplicated_genomes.fa.gz",
        gtdbtk_summary=GTDBTK / "gtdbtk.summary.tsv",
        dram_db=features["databases"]["dram"],
    output:
        annotation=DRAM / "annotate" / "annotations.tsv",
        trnas=DRAM / "annotate" / "trnas.tsv",
        rrnas=DRAM / "annotate" / "rrnas.tsv",
    log:
        DRAM / "annotate.log",
    conda:
        "__environment__.yml"
    container:
        docker["dram"]
    params:
        config=config["dram-config"],
        min_contig_size=1500,
        out_dir=DRAM,
        tmp_dir=DRAM / "annotate",
        parallel_retries=5,
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["quitehighmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["verylongrun"],
        nvme = config["resources"]["nvme"]["small"],
        partition = config["resources"]["partition"]["longrun"]
    shell:
        """
        rm -rf {params.tmp_dir}
        
        echo "Hostname: $(hostname)" 2>> {log} 1>&2
        echo "Temporary directory: $TMPDIR" 2>> {log} 1>&2
        df -h 2>> {log} 1>&2
        
        DRAM.py annotate \
                --config_loc {params.config} \
                --input_fasta {input.dereplicated_genomes} \
                --output_dir {params.tmp_dir} \
                --threads {threads} \
                --gtdb_taxonomy {input.gtdbtk_summary} \
        2>> {log} 1>&2
    """

rule annotate__fix_dram_annotations_scaffold:
    input:
        DRAM / "annotate" / "annotations.tsv",
    output:
        DRAM / "annotate" / "annotations.fixed.tsv",
    log:
        DRAM / "fix_annotation.log",
    container:
        docker["assemble"]
    params:
        script_folder=SCRIPT_FOLDER,
    shell:
        """
        python {params.script_folder}/fix_annotations.py {input} {output} 2>> {log} 1>&2
        """

rule annotate__dram__distill:
    """Distill DRAM annotations."""
    input:
        annotations=DRAM / "annotate" / "annotations.fixed.tsv",
        trnas=DRAM / "annotate" / "trnas.tsv",
        rrnas=DRAM / "annotate" / "rrnas.tsv",
        dram_db=features["databases"]["dram"],
    output:
        genome=DRAM / "genome_stats.tsv",
        metabolism=DRAM / "metabolism_summary.xlsx",
        product_html=DRAM / "product.html",
        product_tsv=DRAM / "product.tsv",
    log:
        DRAM / "distill.log2",
    conda:
        "__environment__.yml"
    container:
        docker["dram"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
    params:
        config=config["dram-config"],
        outdir_tmp=DRAM / "distill",
        outdir=DRAM,
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

rule annotate__dram:
    """Run DRAM on dereplicated genomes."""
    input:
        rules.annotate__dram__distill.output,
