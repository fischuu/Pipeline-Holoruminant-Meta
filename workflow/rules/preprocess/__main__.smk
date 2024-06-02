include: "__functions__.smk"
include: "bowtie2.smk"
include: "fastp.smk"
include: "fastqc.smk"
include: "kraken2.smk"
include: "humann.smk"
include: "metaphlan.smk"
include: "nonpareil.smk"
include: "samtools.smk"
include: "singlem.smk"

rule preprocess:
    """Run the preprocessing steps, included he evaluation ones"""
    input:
        rules.preprocess__bowtie2.input,
        rules.preprocess__fastp.input,
        rules.preprocess__kraken2.input,
        rules.preprocess__humann.input,
        rules.preprocess__metaphlan.input,
        rules.preprocess__nonpareil.input,
        rules.preprocess__samtools.input,
        rules.preprocess__singlem.input,
        