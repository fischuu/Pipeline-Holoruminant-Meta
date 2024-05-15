include: "functions.smk"


# rule vamb_concatenate_one:
#     input:
#         assembly=MEGAHIT_RENAMING / "{assembly_id}.fa",
#     output:
#         concatenated=VAMB / "concatenated" / "{assembly_id}.fa.gz",
#     log:
#         VAMB / "concatenated/{assembly_id}.log",
#     conda:
#         "vamb.yml"
#     shell:
#         """
#         concatenate.py \
#             {output.concatenated} \
#             {input.assembly} \
#         2> {log} 1>&2
#         """


rule vamb_concatenate_one:
    input:
        assembly=ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        concatenated=VAMB / "concatenated" / "{assembly_id}.fa.gz",
    log:
        VAMB / "concatenated/{assembly_id}.log",
    conda:
        "vamb.yml"
    params:
        min_length=params["bin"]["vamb"]["min_length"],
    shell:
        """
        (seqtk seq \
            -L 2000 \
            {input.assembly} \
        | pigz \
        > {output.concatenated} \
        ) 2> {log}
        """


rule vamb_index_one:
    input:
        concatenated=VAMB / "concatenated" / "{assembly_id}.fa.gz",
    output:
        index=touch(VAMB / "indexes" / "{assembly_id}"),
    log:
        VAMB / "indexes/{assembly_id}.log",
    conda:
        "vamb.yml"
    threads: 24
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.concatenated} \
            {output.index} \
        2> {log} 1>&2
        """


rule vamb_map_one:
    input:
        index=VAMB / "indexes" / "{assembly_id}",
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
    output:
        bam=VAMB / "bams" / "{assembly_id}.{sample_id}.{library_id}.bam",
    log:
        VAMB / "bams/{assembly_id}.{sample_id}.{library_id}.log",
    conda:
        "vamb.yml"
    threads: 24
    shell:
        """
        (bowtie2 \
            --threads {threads} \
            --no-unal \
            -x {input.index} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
        | samtools view \
            -F 3584 \
            -u \
        | samtools sort \
            -@ {threads} \
            -l 1 \
            -o {output.bam} \
        ) 2> {log} 1>&2
        """


rule vamb_one:
    input:
        concatenated=VAMB / "concatenated" / "{assembly_id}.fa.gz",
        bams=get_vamb_bams_from_assembly_id,
    output:
        folder=directory(VAMB / "bins" / "{assembly_id}"),
    log:
        VAMB / "bins/{assembly_id}.log",
    conda:
        "vamb.yml"
    params:
        extra="",
    threads: 1
    shell:
        """
        vamb \
            --outdir {output.folder} \
            --fasta {input.concatenated} \
            --bamfiles {input.bams} \
            {params.extra} \
            -p {threads} \
        2> {log} 1>&2
        """


rule vamb:
    input:
        [VAMB / "bins" / f"{assembly_id}" for assembly_id in ASSEMBLIES],
