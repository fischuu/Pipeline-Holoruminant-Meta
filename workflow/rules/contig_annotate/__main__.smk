include: "prodigal.smk"
include: "eggnog.smk"

rule contig_annotate:
    """Annotate on contig level"""
    input:
        expand(rules.contig_assemble__prodigal.input[0], assembly_id=ASSEMBLIES),
        expand(rules.contig_assemble__eggnog.input[0], assembly_id=ASSEMBLIES)
