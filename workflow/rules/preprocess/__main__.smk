include: "__functions__.smk"
include: "bowtie2.smk"
include: "fastp.smk"
include: "fastqc.smk"
include: "samtools.smk"

rule preprocess:
    """Run the preprocessing steps, included he evaluation ones"""
    input:
        rules.preprocess__bowtie2.input,
        rules.preprocess__fastp.input,
        rules.preprocess__samtools.input,
