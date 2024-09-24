def _contig_annotate_aggregate_assembly_eggnog_search(wildcards):
    checkpoint_outputEGG = checkpoints._contigAnnotate__cut_prodigal.get(**wildcards).output[0]
    files = glob_wildcards(os.path.join(checkpoint_outputEGG, "prodigal.chunk.{i}")).i
    return expand(CONTIG_EGGNOG / "{assembly_id}/Chunks/prodigal.chunk.{i}.emapper.seed_orthologs",
                  i=files,
                  assembly_id=wildcards.assembly_id)