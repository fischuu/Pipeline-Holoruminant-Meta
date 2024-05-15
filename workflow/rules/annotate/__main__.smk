include: "quast.smk"
include: "gtdbtk.smk"
include: "dram.smk"
include: "checkm2.smk"


rule annotate:
    """Evaluate the dereplication steps"""
    input:
        rules.annotate__quast.output,
        rules.annotate__checkm2.output,
        rules.annotate__dram.input,
