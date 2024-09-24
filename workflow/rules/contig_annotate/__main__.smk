include: "__functions__.smk"
include: "featurecounts.smk"
include: "eggnog.smk"
include: "prodigal.smk"

rule contig_annotate:
    """Annotate on contig level"""
    input:
        expand(
            rules.contig_assemble__eggnog.input,
            assembly_id=[assembly_id for assembly_id in ASSEMBLIES]
        ),
        expand(
            rules.contig_annotate__featurecount.input,
            assembly_id=[assembly_id for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY],
            sample_id=[sample_id for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY],
            library_id=[library_id for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY]
        ),




