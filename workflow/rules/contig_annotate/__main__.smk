include: "__functions__.smk"
include: "camper.smk"
include: "featurecounts.smk"
include: "eggnog.smk"
include: "prodigal.smk"
include: "sylph.smk"

rule contig_annotate:
    """Annotate on contig level"""
    input:
        rules.contig_annotate__camper.input,
        [
            CONTIG_EGGNOG / f"{assembly_id}/eggnog_output.emapper.annotations"
            for assembly_id in ASSEMBLIES
        ],
        lambda wildcards: [
            CONTIG_FEATURECOUNTS / f"{assembly_id}/{assembly_id}_{sample_id}.{library_id}.prodigal_fc.txt"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
        rules.contig_annotate__sylph.input,
        

