# ---------- preprocess__fastp__run ----------
rule preprocess__fastp__run:
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
        FASTP / "benchmark/{sample_id}.{library_id}.tsv",
    container:
        docker["preprocess"],
    params:
        adapter_forward=get_forward_adapter,
        adapter_reverse=get_reverse_adapter,
        extra=params["preprocess"]["fastp"]["extra"],
        length_required=params["preprocess"]["fastp"]["length_required"],
        temp_forward_=lambda w: FASTP / f"{w.sample_id}.{w.library_id}_tmp_1.fq",
        temp_reverse_=lambda w: FASTP / f"{w.sample_id}.{w.library_id}_tmp_2.fq",
        temp_unpaired1=lambda w: FASTP / f"{w.sample_id}.{w.library_id}_tmp_u1.fq",
        temp_unpaired2=lambda w: FASTP / f"{w.sample_id}.{w.library_id}_tmp_u2.fq",
    threads: esc("cpus", "preprocess__fastp__run"),
    resources:
        runtime=esc("runtime", "preprocess__fastp__run"),
        mem_mb=esc("mem_mb", "preprocess__fastp__run"),
        cpu_per_task=esc("cpus", "preprocess__fastp__run"),
        partition=esc("partition", "preprocess__fastp__run"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "preprocess__fastp__run", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__fastp__run")),
    shell: """
        echo "LOCAL_SCRATCH is " $LOCAL_SCRATCH 2> {log}.{resources.attempt} 1>&2

        fastp \
            --in1 {input.forward_} \
            --in2 {input.reverse_} \
            --out1 {params.temp_forward_} \
            --out2 {params.temp_reverse_} \
            --unpaired1 {params.temp_unpaired1} \
            --unpaired2 {params.temp_unpaired2} \
            --html {output.html} \
            --json {output.json} \
            --verbose \
            --adapter_sequence {params.adapter_forward} \
            --adapter_sequence_r2 {params.adapter_reverse} \
            --length_required {params.length_required} \
            --thread {threads} \
            {params.extra} \
        2> {log}.{resources.attempt} 1>&2

        bgzip -l 9 -@ {threads} {params.temp_forward_}
        bgzip -l 9 -@ {threads} {params.temp_reverse_}
        bgzip -l 9 -@ {threads} {params.temp_unpaired1}
        bgzip -l 9 -@ {threads} {params.temp_unpaired2}

        mv {params.temp_forward_}.gz {output.forward_}
        mv {params.temp_reverse_}.gz {output.reverse_}
        mv {params.temp_unpaired1}.gz {output.unpaired1}
        mv {params.temp_unpaired2}.gz {output.unpaired2}

        echo "Checking integrity of output files" >> {log}.{resources.attempt}
        gzip -t {output.forward_} 2>> {log}.{resources.attempt}
        gzip -t {output.reverse_} 2>> {log}.{resources.attempt}
        gzip -t {output.unpaired1} 2>> {log}.{resources.attempt}
        gzip -t {output.unpaired2} 2>> {log}.{resources.attempt}
        echo "Integrity check completed" >> {log}.{resources.attempt}

        mv {log}.{resources.attempt} {log}
    """

# ---------- preprocess__fastp ----------
rule preprocess__fastp:
    """Aggregate all fastp outputs"""
    input:
        reads=[
            FASTP / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in "1 2 u1 u2".split()
        ],
        json=[
            FASTP / f"{sample_id}.{library_id}_fastp.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
