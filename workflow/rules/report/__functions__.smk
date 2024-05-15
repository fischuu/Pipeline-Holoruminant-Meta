def _path_to_string(path_list):
    """Convert a list of paths to a list of strings"""
    return [str(x) for x in path_list]


def get_stats_files_from_sample_and_library_ids(wildcards):
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id

    reads_fastqc = [
        READS / f"{sample_id}.{library_id}_{end}_fastqc.zip" for end in ["1", "2"]
    ]

    pre_fastp_fastqc = [
        FASTP / f"{sample_id}.{library_id}_{end}_fastqc.zip"
        for end in ["1", "2", "u1", "u2"]
    ]

    pre_fastp_html = FASTP / f"{sample_id}.{library_id}_fastp.html"

    pre_bowtie2 = (
        [
            PRE_BOWTIE2 / genome / f"{sample_id}.{library_id}.{extension}"
            for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
            for genome in HOST_NAMES
        ]
        if len(HOST_NAMES) > 0
        else []
    )

    pre_nonhost_fastqc = (
        [
            PRE_BOWTIE2 / f"non{genome}" / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for end in ["1", "2"]
            for genome in HOST_NAMES
        ]
        if len(HOST_NAMES) > 0
        else []
    )

    pre_kraken2 = [
        KRAKEN2 / kraken_db / f"{sample_id}.{library_id}.report"
        for kraken_db in KRAKEN2_DBS
    ]

    quantify_bowtie2 = [
        QUANT_BOWTIE2 / f"{sample_id}.{library_id}.{extension}"
        for extension in ["stats.txt", "flagstats.txt"]
    ]

    all_files = (
        _path_to_string(reads_fastqc)
        + _path_to_string(pre_fastp_fastqc)
        + [pre_fastp_html]
        + _path_to_string(pre_bowtie2)
        + _path_to_string(pre_nonhost_fastqc)
        + _path_to_string(pre_kraken2)
        + _path_to_string(quantify_bowtie2)
    )

    return all_files
