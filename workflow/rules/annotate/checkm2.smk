rule _annotate__checkm2__predict:
    """Run CheckM2 over the dereplicated mags"""
    input:
        mags=DREP / "dereplicated_genomes",
        db=features["databases"]["checkm2"],
    output:
        CHECKM / "quality_report.tsv",
    log:
        CHECKM / "quality_report.log",
    threads: 24
    conda:
        "checkm2.yml"
    singularity:
        docker["checkm2"]
    params:
        out_dir=CHECKM / "predict",
    resources:
        mem_mb=16 * 1024,
    shell:
        """
        rm -rfv {params.out_dir} 2> {log} 1>&2

        checkm2 predict \
            --threads {threads} \
            --input {input.mags} \
            --extension .fa.gz \
            --output-directory {params.out_dir} \
            --database_path {input.db}/uniref100.KO.1.dmnd \
            --remove_intermediates \
        2>> {log} 1>&2

        mv {params.out_dir}/quality_report.tsv {output} 2>> {log} 1>&2
        rm --recursive --verbose --force {params.out_dir} 2>> {log} 1>&2
        """


rule annotate__checkm2:
    """Run CheckM2"""
    input:
        rules._annotate__checkm2__predict.output,
