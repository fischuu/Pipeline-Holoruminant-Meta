include: "__functions__.smk"
include: "bowtie2.smk"
include: "diamond.smk"
include: "fastp.smk"
include: "fastqc.smk"
include: "kraken2.smk"
include: "krona.smk"
include: "humann.smk"
include: "metaphlan.smk"
include: "phyloflash.smk"
include: "nonpareil.smk"
include: "samtools.smk"
include: "singlem.smk"
include: "sylph.smk"

rule preprocess:
    """Run the preprocessing steps, included he evaluation ones"""
    input:
        rules.preprocess__bowtie2.input,
        rules.preprocess__diamond.input,
        rules.preprocess__fastp.input,
        rules.preprocess__kraken2.input,
        rules.preprocess__krona.input,
        rules.preprocess__humann.input,
        rules.preprocess__metaphlan.input,
        rules.preprocess__nonpareil.input,
        rules.preprocess__phyloflash.input,
        rules.preprocess__samtools.input,
        rules.preprocess__singlem.input,
        rules.preprocess__sylph.input,
