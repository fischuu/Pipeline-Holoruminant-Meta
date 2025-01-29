rule _annotate__camper__annotate:
    """Annotate dereplicate genomes with CAMPER"""
    input:
        dereplicated_genomes=DREP / "dereplicated_genomes.fa.gz"
    output:
        annotation=CAMPER / "annotations.tsv"
    log:
        CAMPER / "annotate.log",
    conda:
        "__environment__.yml"
    container:
        docker["camper"]
    params:
        out_dir=CAMPER,
        parallel_retries=5,
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
        nvme = config["resources"]["nvme"]["small"]
    shell:
        """
        camper_annotate -i {input} \
                        -o {output} \
                        --threads {threads} \

        camper_distill  -a <path to annotations.tsv> -o <name of output.tsv>

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

rule annotate__camper:
    """Run CAMPER on dereplicated genomes."""
    input:
        rules._annotate__camper__annotate,
