rule _preprocess__humann__run:
    """Run HumanN3 over one sample
    """
    input:
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        prot_dbs=features["databases"]["humann_prot_dbs"],
        nt_dbs=features["databases"]["humann_nt_dbs"],
    output:
        humann_out=HUMANN / "{sample_id}.{library_id}_genefamilies.tsv",
    log:
        HUMANN / "log" / "{sample_id}.{library_id}.log",
    benchmark:
        HUMANN / "benchmark" / "{sample_id}.{library_id}.tsv",
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
    params:
        out_folder=HUMANN,
        out_name="{sample_id}.{library_id}",
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"]
    shell:
        """
        cat {input.forward_} {input.reverse_} | \
        humann --input - --output {params.out_folder} \
        --protein-database {input.prot_dbs} \
        --nucleotide-database {input.nt_dbs} \
        --output-basename {params.out_name}

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
    singularity:
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