rule _reads__link:
    """Make a link to the original file, with a prettier name than default"""
    input:
        forward_=get_forward,
        reverse_=get_reverse,
    output:
        forward_=READS / "{sample}.{library}_1.fq.gz",
        reverse_=READS / "{sample}.{library}_2.fq.gz",
    log:
        READS / "{sample}.{library}.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["reads"]
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_} 2>  {log} 1>&2
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_} 2>> {log} 1>&2
        """


rule reads__link:
    """Link all reads in the samples.tsv"""
    input:
        [
            READS / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],
