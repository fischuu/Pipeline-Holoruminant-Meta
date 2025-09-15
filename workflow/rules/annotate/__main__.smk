include: "bakta.smk"
include: "quast.smk"
include: "gtdbtk.smk"
include: "dram.smk"
include: "dram_mag.smk"
include: "eggnog.smk"
include: "checkm2.smk"
include: "proteinortho.smk"
include: "phylophlan.smk"


rule annotate:
    """Evaluate the dereplication steps"""
    input:
        rules.annotate__bakta.output,
        rules.annotate__bakta_mags.input,
        rules.annotate__eggnog.output,
        rules.annotate__quast.output,
        rules.annotate__checkm2__predict.output,
        rules.annotate__gtdbtk__classify.output,
        rules.annotate__dram.input,
        rules.annotate__dram_mags.input,
        rules.annotate__proteinortho.output,
        rules.annotate__phylophlan.output,
