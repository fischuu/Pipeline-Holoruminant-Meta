rule _annotate__eggnog:
    """Run eggnog over the dereplicated mags"""
    input:
        contigs=DREP / "dereplicated_genomes.fa.gz",
    output:
        directory(EGGNOG) ,
    log:
        EGGNOG / "eggnog.log",
    singularity:
        docker["eggnog"]
    params:
        out_dir=EGGNOG,
        db=features["databases"]["eggnog"],
        prefix="eggnog"
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"]//config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
        nvme = config["resources"]["nvme"]["large"]
    shell:
        """
        cp {params.db}/* $TMPDIR  2>> {log} 1>&2;
        
        emapper.py -m diamond \
                   --data_dir $TMPDIR \
                   --itype metagenome \
                   --genepred prodigal \
                   --dbmem \
                   --no_annot \
                   --no_file_comments \
                   --cpu {threads} \
                   -i {input.contigs} \
                   --output_dir {params.out_dir} \
                   -o {params.prefix}  \
                   2>> {log} 1>&2;
        """
