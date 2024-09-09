rule _assemble__metaspades:
    """Run MetaSPAdes over one sample, merging all libraries in the process"""
    input:
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        fasta=METASPADES / "{assembly_id}.fa.gz",
        tarball=METASPADES / "{assembly_id}.tar.gz",
    log:
        log=METASPADES / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        attempt=get_attempt,
    params:
        out_dir=lambda w: METASPADES / w.assembly_id,
        kmer_size=params["assemble"]["metaspades"]["kmer_size"],
        min_contig_len=params["assemble"]["metaspades"]["min_contig_len"],
        forwards=aggregate_forwards_for_metaspades,
        reverses=aggregate_reverses_for_metaspades,
        assembly_id=lambda w: w.assembly_id,
    shell:
        """
        metaspades.py \
            -t {threads} \
            -m {resources.mem_per_cpu * threads} \
            -k {params.kmer_size} \
            --only-assembler \
            -1 {params.forwards} \
            -2 {params.reverses} \
            -o {params.out_dir} \
        2> {log}.{resources.attempt} 1>&2

        ( seqtk seq \
            {params.out_dir}/contigs.fasta \
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
