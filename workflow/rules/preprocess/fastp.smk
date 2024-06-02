rule _preprocess__fastp__run:
    """Run fastp on one library"""
    input:
        forward_=READS / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=READS / "{sample_id}.{library_id}_2.fq.gz",
    output:
        forward_=FASTP / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=FASTP / "{sample_id}.{library_id}_2.fq.gz",
        unpaired1=FASTP / "{sample_id}.{library_id}_u1.fq.gz",
        unpaired2=FASTP / "{sample_id}.{library_id}_u2.fq.gz",
        html=FASTP / "{sample_id}.{library_id}_fastp.html",
        json=FASTP / "{sample_id}.{library_id}_fastp.json",
    log:
        FASTP / "{sample_id}.{library_id}.log",
    benchmark:
        FASTP / "benchmark/{sample_id}.{library_id}.tsv"
    conda:
        "__environment__.yml"
    singularity:
        docker["preprocess"]
    params:
        adapter_forward=get_forward_adapter,
        adapter_reverse=get_reverse_adapter,
        extra=params["preprocess"]["fastp"]["extra"],
        length_required=params["preprocess"]["fastp"]["length_required"],
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"]
    shell:
        """
        fastp \
            --in1 {input.forward_} \
            --in2 {input.reverse_} \
            --out1 >(bgzip -l 9 -@ {threads} > {output.forward_}) \
            --out2 >(bgzip -l 9 -@ {threads} > {output.reverse_}) \
            --unpaired1 >(bgzip -l 9 -@ {threads} > {output.unpaired1}) \
            --unpaired2 >(bgzip -l 9 -@ {threads} > {output.unpaired2}) \
            --html {output.html} \
            --json {output.json} \
            --verbose \
            --adapter_sequence {params.adapter_forward} \
            --adapter_sequence_r2 {params.adapter_reverse} \
            --length_required {params.length_required} \
            --thread {threads} \
            {params.extra} \
        2> {log} 1>&2
        """


rule preprocess__fastp:
    """Get all files from fastp"""
    input:
        reads=[
            FASTP / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in "1 2 u1 u2".split(" ")
        ],
        json=[
            FASTP / f"{sample_id}.{library_id}_fastp.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
