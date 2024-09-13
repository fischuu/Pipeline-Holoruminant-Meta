include: "prodigal.smk"

rule contig_annotate:
    """Annotate on contig level"""
    input:
        expand(rules.contig_assemble__prodigal.input[0], assembly_id=ASSEMBLIES)
