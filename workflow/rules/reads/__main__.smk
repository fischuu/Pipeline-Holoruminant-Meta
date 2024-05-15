include: "__functions__.smk"
include: "link.smk"
include: "fastqc.smk"


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads__link.input,
        rules.reads__fastqc.input,
