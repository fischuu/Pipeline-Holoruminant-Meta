include: "__functions__.smk"
include: "bowtie2.smk"
include: "concoct.smk"
include: "drep.smk"
include: "magscot.smk"
include: "maxbin2.smk"
include: "megahit.smk"
include: "metabat2.smk"


rule assemble:
    """Run the assemble module"""
    input:
        rules.assemble__drep.input,
