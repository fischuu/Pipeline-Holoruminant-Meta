include: "__functions__.smk"
include: "bowtie2.smk"
include: "coverm.smk"
include: "samtools.smk"


rule quantify:
    input:
        rules.quantify__coverm.input,
        rules.quantify__samtools.input,
