include: "quast.smk"
include: "gtdbtk.smk"
include: "dram.smk"
include: "checkm2.smk"


rule annotate:
    """Evaluate the dereplication steps"""
    input:
        rules.annotate__quast.output,
        rules._annotate__checkm2__predict.output,
        rules._annotate__gtdbtk__classify.output,
 #       rules.annotate__dram.input,
