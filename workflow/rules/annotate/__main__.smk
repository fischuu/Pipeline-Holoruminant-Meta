include: "bakta.smk"
include: "quast.smk"
include: "gtdbtk.smk"
include: "dram.smk"
include: "eggnog.smk"
include: "checkm2.smk"
include: "proteinortho.smk"
include: "phylophlan.smk"

rule annotate:
    """Evaluate the dereplication steps"""
    input:
        rules._annotate__bakta.output,
        rules._annotate__eggnog.output,
        rules.annotate__quast.output,
        rules._annotate__checkm2__predict.output,
        rules._annotate__gtdbtk__classify.output,
        rules.annotate__dram.input,
        rules._annotate__proteinortho_new.output,
        rules._annotate__phylophlan.input,