rule _assemble__megahit:
    """Run megahit over one sample, merging all libraries in the process

    Note: the initial rm -rf is to delete the folder that snakemake creates.
    megahit refuses to overwrite an existing folder
    """
    input:
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        fasta=MEGAHIT / "{assembly_id}.fa.gz",
        tarball=MEGAHIT / "{assembly_id}.tar.gz",
    log:
        log=MEGAHIT / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
        attempt=get_attempt,
    params:
        out_dir=lambda w: MEGAHIT / w.assembly_id,
        mincount=params["assemble"]["megahit"]["mincount"],
        kmin=params["assemble"]["megahit"]["kmin"],
        kmax=params["assemble"]["megahit"]["kmax"],
        kstep=params["assemble"]["megahit"]["kstep"],
        min_contig_len=params["assemble"]["megahit"]["min_contig_len"],
        forwards=aggregate_forwards_for_megahit,
        reverses=aggregate_reverses_for_megahit,
        assembly_id=lambda w: w.assembly_id,
    shell:
        """
        megahit \
            --num-cpu-threads {threads} \
            --min-count {params.mincount} \
            --k-min {params.kmin} \
            --k-max {params.kmax} \
            --k-step {params.kstep} \
            --min-contig-len {params.min_contig_len} \
            --verbose \
            --force \
            --out-dir {params.out_dir} \
            --continue \
            -1 {params.forwards} \
            -2 {params.reverses} \
        2> {log}.{resources.attempt} 1>&2

        ( seqtk seq \
            {params.out_dir}/final.contigs.fa \
        | cut -f 1 -d " " \
        | paste - - \
        | awk \
            '{{printf(">{params.assembly_id}:bin_NA@contig_%08d\\n%s\\n", NR, $2)}}' \
        | bgzip \
            -l9 \
            -@ {threads} \
        > {output.fasta} \
        ) 2>> {log}.{resources.attempt}

        tar \
            --create \
            --file {output.tarball} \
            --remove-files \
            --use-compress-program="pigz --best --processes {threads}" \
            --verbose \
            {params.out_dir} \
        2>> {log}.{resources.attempt} 1>&2

        mv {log}.{resources.attempt} {log}
        """


rule assemble__megahit:
    """Rename all assemblies contigs to avoid future collisions"""
    input:
        [MEGAHIT / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
