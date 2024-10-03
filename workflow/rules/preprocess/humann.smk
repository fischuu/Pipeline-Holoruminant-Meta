rule _preprocess__humann__run:
    """Run HumanN3 over one sample
    """
    input:
        forward_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        prot_dbs=features["databases"]["humann_prot_dbs"],
        nt_dbs=features["databases"]["humann_nt_dbs"],
        mp_out=METAPHLAN / "profiled" / "{sample_id}.{library_id}.txt",
    output:
        humann_out=HUMANN / "{sample_id}.{library_id}_genefamilies.tsv",
        cat= temp(HUMANN / "{sample_id}.{library_id}_concatenated.fastq.gz"),
    log:
        HUMANN / "log" / "{sample_id}.{library_id}.log",
    benchmark:
        HUMANN / "benchmark" / "{sample_id}.{library_id}.tsv",
    conda:
        "__environment__.yml"
    container:
        docker["preprocess"]
    params:
        out_folder=HUMANN,
        out_name="{sample_id}.{library_id}",
        additional_options=params["preprocess"]["humann"]["additional_options"],
        tmp = config["tmp_storage"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"]
    shell:
        """
        TMPDIR={params.tmp}

        cat {input.forward_} {input.reverse_} > {output.cat}
        
        humann --input {output.cat} --output {params.out_folder} \
        --threads {threads} \
        --protein-database {input.prot_dbs} \
        --nucleotide-database {input.nt_dbs} \
        --output-basename {params.out_name} \
        --taxonomic-profile {input.mp_out} \
        {params.additional_options} \
        2>> {log} 1>&2

        """


rule _preprocess__humann__condense:
    """Aggregate all the HumanN results into a single table"""
    input:
        genefamily_data=[
            HUMANN / f"{sample_id}.{library_id}_genefamilies.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
    output:
        HUMANN / "humann_genefamilies.tsv",
    log:
        HUMANN / "humann.log",
    benchmark:
        HUMANN / "benchmark/humann.tsv",
    conda:
        "__environment__.yml"
    container:
        docker["preprocess"]
    params:
        input_dir=HUMANN,
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"]
    shell:
        """
        cat {input.genefamily_data} > {output}
        2> {log} 1>&2
        """

rule preprocess__humann:
    input:
        rules._preprocess__humann__condense.output