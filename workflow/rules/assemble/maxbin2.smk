rule _assemble__maxbin2__run:
    """Run MaxBin2 over a single assembly"""
    input:
        assembly=MEGAHIT / "{assembly_id}.fa.gz",
        crams=get_crams_from_assembly_id,
    output:
        workdir=directory(MAXBIN2 / "{assembly_id}"),
    log:
        MAXBIN2 / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
    params:
        seed=1,
        coverage=lambda w: MAXBIN2 / f"{w.assembly_id}/maxbin2.coverage",
    shell:
        """
        mkdir --parents {output.workdir}

        ( samtools coverage {input.crams} \
        | awk '{{print $1"\\t"$5}}' \
        | grep -v '^#' \
        ) > {params.coverage} \
        2> {log}

        run_MaxBin.pl \
            -thread {threads} \
            -contig {input.assembly} \
            -out {output.workdir}/maxbin2 \
            -abund {params.coverage} \
        2> {log} 1>&2

        rename \
            's/\\.fasta$/.fa/' \
            {output.workdir}/*.fasta \
        2>> {log}

        find \
            {output.workdir} \
            -name "*.fa" \
            -exec pigz --best --verbose {{}} \; \
        2>> {log} 1>&2

        rm \
            --recursive \
            --force \
            {output.workdir}/maxbin.{{coverage,log,marker,noclass,summary,tooshort}} \
            {output.workdir}/maxbin2.marker_of_each_bin.tar.gz \
        2>> {log} 1>&2
        """


rule assemble__maxbin2:
    """Run MaxBin2 over all assemblies"""
    input:
        [MAXBIN2 / assembly_id for assembly_id in ASSEMBLIES],
