rule preprocess__fastqc__fastp:
    """Run fastqc over all libraries after fastp"""
    input:
        [
            FASTP / f"{sample_id}.{library_id}_{end}_fastqc.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in "1 2".split(" ")
            for extension in "html zip".split(" ")
        ],


rule preprocess__fastqc__nonhost:
    """Run fastqc over all libraries after fastp"""
    input:
        [
            PRE_BOWTIE2
            / f"non{genome}/{sample_id}.{library_id}_{end}_fastqc.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in "1 2".split(" ")
            for extension in "html zip".split(" ")
            for genome in HOST_NAMES
        ],


rule preprocess__fastqc:
    """Run fastqc over all fastq files in prepreprocess"""
    input:
        rules.preprocess__fastqc__fastp.input,
        rules.preprocess__fastqc__nonhost.input,
