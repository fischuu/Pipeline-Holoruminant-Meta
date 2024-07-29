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
    shell:        """
        # Define intermediate and final output file paths
        intermediate_out1=$(mktemp {output.forward_}.XXXXXX)
        intermediate_out2=$(mktemp {output.reverse_}.XXXXXX)
        intermediate_unpaired1=$(mktemp {output.unpaired1}.XXXXXX)
        intermediate_unpaired2=$(mktemp {output.unpaired2}.XXXXXX)

        # Run fastp with intermediate files
        fastp \
            --in1 {input.forward_} \
            --in2 {input.reverse_} \
            --out1 "$intermediate_out1" \
            --out2 "$intermediate_out2" \
            --unpaired1 "$intermediate_unpaired1" \
            --unpaired2 "$intermediate_unpaired2" \
            --html {output.html} \
            --json {output.json} \
            --verbose \
            --adapter_sequence {params.adapter_forward} \
            --adapter_sequence_r2 {params.adapter_reverse} \
            --length_required {params.length_required} \
            --thread {threads} \
            {params.extra} \
        2> {log} 1>&2

        # Compress the outputs using bgzip
        bgzip -l 9 -@ {threads} "$intermediate_out1"
        bgzip -l 9 -@ {threads} "$intermediate_out2"
        bgzip -l 9 -@ {threads} "$intermediate_unpaired1"
        bgzip -l 9 -@ {threads} "$intermediate_unpaired2"

        # Move the compressed files to the final destination
        mv "$intermediate_out1.gz" {output.forward_}
        mv "$intermediate_out2.gz" {output.reverse_}
        mv "$intermediate_unpaired1.gz" {output.unpaired1}
        mv "$intermediate_unpaired2.gz" {output.unpaired2}

        # Check the integrity of the gzipped files and log the output
        echo "Checking integrity of output files" >> {log}
        gzip -t {output.forward_} 2>> {log}
        gzip -t {output.reverse_} 2>> {log}
        gzip -t {output.unpaired1} 2>> {log}
        gzip -t {output.unpaired2} 2>> {log}
        echo "Integrity check completed" >> {log}
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
