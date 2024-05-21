rule _assemble__drep__separate_bins:
    input:
        assemblies=[MAGSCOT / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
    output:
        out_dir=directory(DREP / "separated_bins"),
    log:
        DREP / "separate_bins.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["assemble"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        mkdir --parents {output.out_dir} 2> {log} 1>&2

        ( gzip \
            --decompress \
            --stdout \
            {input.assemblies} \
        | paste - - \
        | tr -d ">" \
        | tr "@" "\t" \
        | awk \
            '{{print ">" $1 "@" $2 "\\n" $3 > "{output.out_dir}/" $1 ".fa" }}' \
        ) >> {log} 2>&1
        """


rule _assemble__drep__run:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=DREP / "separated_bins",
    output:
        dereplicated_genomes=directory(DREP / "dereplicated_genomes"),
        data=DREP / "data.tar.gz",
        data_tables=DREP / "data_tables.tar.gz",
    log:
        DREP / "drep.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
        attempt=get_attempt,
    params:
        out_dir=DREP,
    shell:
        """
        rm \
            --recursive \
            --force \
            {params.out_dir}/data_tables \
            {params.out_dir}/data \
            {params.out_dir}/dereplicated_genomes \
            {params.out_dir}/figures \
            {params.out_dir}/log \
        2> {log}.{resources.attempt} 1>&2

        dRep dereplicate \
            {params.out_dir} \
            --processors {threads} \
            --completeness 50 \
            --S_ani 0.9 \
            --genomes {input.genomes}/*.fa \
        2>> {log}.{resources.attempt} 1>&2

        for folder in data data_tables ; do
            tar \
                --create \
                --directory {params.out_dir} \
                --file {params.out_dir}/${{folder}}.tar.gz \
                --remove-files \
                --use-compress-program="pigz --processes {threads}" \
                --verbose \
                ${{folder}} \
            2>> {log} 1>&2
        done

        mv {log}.{resources.attempt} {log}
        """


rule _assemble__drep__join_genomes:
    """Join all the dereplicated genomes into a single file."""
    input:
        DREP / "dereplicated_genomes",
    output:
        DREP / "dereplicated_genomes.fa.gz",
    log:
        DREP / "dereplicated_genomes.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        ( zcat \
            {input}/*.fa.gz \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output} \
        ) 2> {log}
        """


rule assemble__drep:
    input:
        DREP / "dereplicated_genomes.fa.gz",
