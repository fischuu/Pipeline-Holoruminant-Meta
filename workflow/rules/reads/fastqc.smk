rule reads__fastqc:
    """Get all fastqc reports of the raw reads"""
    input:
        [
            READS / f"{sample}.{library}_{end}_fastqc.zip"
            for sample, library in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],
