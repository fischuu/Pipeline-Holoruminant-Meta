rule _annotate__dram__annotate:
    """Annotate dereplicate genomes with DRAM"""
    input:
        dereplicated_genomes=DREP / "dereplicated_genomes.fa.gz",
        gtdbtk_summary=GTDBTK / "gtdbtk.summary.tsv",
        dram_db=features["databases"]["dram"],
    output:
        annotation=DRAM / "annotations.tsv",
        trnas=DRAM / "trnas.tsv",
        rrnas=DRAM / "rrnas.tsv",
        tarball=DRAM / "annotate.tar.gz",
    log:
        DRAM / "annotate.log",
    conda:
        "__environment__.yml"
    singularity:
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
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
#        nvme = config["resources"]["nvme"]["small"]
    shell:
        """
        rm -rf {params.tmp_dir}
        
        DRAM.py annotate \
                --config_loc {params.config} \
                --input_fasta {input.dereplicated_genomes} \
                --output_dir {params.tmp_dir} \
                --threads {threads} \
                --gtdb_taxonomy {input.gtdbtk_summary} \
        2>> {log} 1>&2

        for file in annotations trnas rrnas ; do
            ( csvstack \
                --tabs \
                {params.tmp_dir}/*/$file.tsv \
            | csvformat \
                --out-tabs \
            > {params.out_dir}/$file.tsv \
            ) 2>> {log}
        done

        tar \
            --create \
            --directory {params.out_dir} \
            --file {output.tarball} \
            --remove-files \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            annotate \
        2>> {log} 1>&2
        """


rule _annotate__dram__distill:
    """Distill DRAM annotations."""
    input:
        annotations=DRAM / "annotations.tsv",
        trnas=DRAM / "trnas.tsv",
        rrnas=DRAM / "rrnas.tsv",
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
    singularity:
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
        rules._annotate__dram__distill.output,
