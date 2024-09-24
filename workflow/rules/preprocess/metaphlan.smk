rule _preprocess__metaphlan__run:
    """Run metaphlan over one sample
    """
    input:
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        data=features["databases"]["metaphlan4"],
    output:
        bt2_out=METAPHLAN / "bowtie2_out" / "{sample_id}.{library_id}.bz2",
        mp_out=METAPHLAN / "profiled" / "{sample_id}.{library_id}.txt",
    log:
        METAPHLAN / "log" / "{sample_id}.{library_id}.log",
    benchmark:
        METAPHLAN / "benchmark" / "{sample_id}.{library_id}.tsv",
    conda:
        "__environment__.yml"
    container:
        docker["preprocess"]
    threads: config["resources"]["cpu_per_task"]["single_thread"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"]
    shell:
        """
        metaphlan {input.forward_},{input.reverse_} \
                  --bowtie2out {output.bt2_out} \
                  --nproc {threads} \
                  --input_type fastq \
                  -o {output.mp_out} \
                  --bowtie2db {input.data} \
        2>> {log} 1>&2
        """


rule _preprocess__metaphlan__condense:
    """Aggregate all the metaphlan results into a single table"""
    input:
        profiled_data=[
            METAPHLAN / "profiled" / f"{sample_id}.{library_id}.txt"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
    output:
        METAPHLAN / "metaphlan_profiled.tsv",
    log:
        METAPHLAN / "metaphlan.log",
    benchmark:
        METAPHLAN / "benchmark/metaphlan.tsv",
    conda:
        "__environment__.yml"
    container:
        docker["preprocess"]
    params:
        input_dir=METAPHLAN,
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"]
    shell:
        """
         merge_metaphlan_tables.py {input.profiled_data} > {output}
        #cat {input.profiled_data} > {output}
        #2> {log} 1>&2
        """

rule preprocess__metaphlan:
    input:
        rules._preprocess__metaphlan__condense.output
