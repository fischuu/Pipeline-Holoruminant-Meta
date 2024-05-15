def get_vamb_bams_from_assembly_id(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    bam_files = []
    for sample_id, library_id in samples_in_assembly:
        bam_files.append(
            VAMB / "bams" / f"{assembly_id}.{sample_id}.{library_id}.bam",
        )
    return bam_files
