rule assemble__maxbin2__run:
    """Run MaxBin2 over a single assembly"""
    input:
        assembly=lambda wildcards: (
            MEGAHIT / f"{wildcards.assembly_id}.fa.gz" if config["assembler"] == "megahit" else
            METASPADES / f"{wildcards.assembly_id}.fa.gz"
        ),
        crams=get_crams_from_assembly_id,
    output:
        workdir=directory(MAXBIN2 / "{assembly_id}"),
    log:
        MAXBIN2 / "{assembly_id}.log",
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__maxbin2__run")
    resources:
        runtime=esc("runtime", "assemble__maxbin2__run"),
        mem_mb=esc("mem_mb", "assemble__maxbin2__run"),
        cpus_per_task=esc("cpus", "assemble__maxbin2__run"),
        slurm_partition=esc("partition", "assemble__maxbin2__run"),
        slurm_extra="'--gres=nvme:" + str(esc_val("nvme", "assemble__maxbin2__run", attempt=1)) + "'",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__maxbin2__run"))
    params:
        seed=1,
        coverage=lambda w: MAXBIN2 / f"{w.assembly_id}/maxbin2.coverage",
        minLen=params["assemble"]["maxbin"]["min_contig_len"],
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
            -min_contig_length {params.minLen} \
        2> {log} 1>&2

        rename \
            's/\\.fasta$/.fa/' \
            {output.workdir}/*.fasta \
        2>> {log}

        fa_files=$(find {output.workdir} -name "*.fa")
        for fa in $fa_files; do
            pigz --best --verbose "$fa"
        done 2>> {log} 1>&2

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
