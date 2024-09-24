rule _assemble__metaspades:
    """Run MetaSPAdes over one sample, merging all libraries in the process"""
    input:
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        fasta=METASPADES / "{assembly_id}.fa.gz",
        tarball=METASPADES / "{assembly_id}.tar.gz",
        concatenated_forwards=temp(METASPADES / "{assembly_id}_R1_concat.fastq.gz"),
        concatenated_reverses=temp(METASPADES / "{assembly_id}_R2_concat.fastq.gz"),
    log:
        log=METASPADES / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    container:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["veryhighmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["verylongrun"],
        partition=config["resources"]["partition"]["highlong"],
        attempt=get_attempt,
    params:
        out_dir=lambda w: METASPADES / w.assembly_id,
        kmer_size=params["assemble"]["metaspades"]["kmer_size"],
        additional_options=params["assemble"]["metaspades"]["additional_options"],
        forwards=aggregate_forwards_for_metaspades,
        reverses=aggregate_reverses_for_metaspades,
        assembly_id=lambda w: w.assembly_id,
    shell:
        """
        # Concatenate forward reads into a single file
        cat {params.forwards} > {output.concatenated_forwards}
        
        # Concatenate reverse reads into a single file
        cat {params.reverses} > {output.concatenated_reverses}
        
        metaspades.py \
            -t {threads} \
            -k {params.kmer_size} \
            {params.additional_options} \
            -1 {output.concatenated_forwards} \
            -2 {output.concatenated_reverses} \
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

rule assemble__metaspades:
    """Rename all assemblies contigs to avoid future collisions"""
    input:
        [METASPADES / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],