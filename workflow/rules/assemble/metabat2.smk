rule assemble__metabat2__run:
    """Run metabat2 end-to-end on a single assembly"""
    input:
        crams=get_crams_from_assembly_id,
        crais=get_crais_from_assembly_id,
        assembly=lambda wildcards: (
            MEGAHIT / f"{wildcards.assembly_id}.fa.gz" if config["assembler"] == "megahit" else
            METASPADES / f"{wildcards.assembly_id}.fa.gz"
        ),
    output:
        bins=directory(METABAT2 / "{assembly_id}"),
    log:
        METABAT2 / "{assembly_id}.log",
    container:
        docker["assemble"]
    params:
        bins_prefix=lambda w: METABAT2 / f"{w.assembly_id}/bin",
        bams=compose_bams_for_metabat2_run,
        depth=lambda w: METABAT2 / f"{w.assembly_id}.depth",
        paired=lambda w: METABAT2 / f"{w.assembly_id}.paired",
        workdir=METABAT2,
        minLen=params["assemble"]["metabat"]["min_contig_len"],
    threads: esc("cpus", "assemble__metabat2__run")
    resources:
        runtime=esc("runtime", "assemble__metabat2__run"),
        mem_mb=esc("mem_mb", "assemble__metabat2__run"),
        cpus_per_task=esc("cpus", "assemble__metabat2__run"),
        slurm_partition=esc("partition", "assemble__metabat2__run"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'assemble__metabat2__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__metabat2__run"))
    shell:
        """
        for cram in {input.crams} ; do

            bam={params.workdir}/$(basename $cram .cram).bam

            samtools view \
                --exclude-flags 4 \
                --fast \
                --output $bam \
                --output-fmt BAM \
                --reference {input.assembly} \
                --threads {threads} \
                $cram

        done 2> {log} 1>&2

        jgi_summarize_bam_contig_depths \
            --outputDepth {params.depth} \
            --pairedContigs {params.paired} \
            {params.bams} \
        2>> {log} 1>&2

        metabat2 \
            --inFile {input.assembly} \
            --abdFile {params.depth} \
            --outFile {params.bins_prefix} \
            --numThreads {threads} \
            --minContig {params.minLen} \
            --verbose \
        2> {log} 1>&2

        rm \
            --force \
            --verbose \
            {params.bams} \
            {params.depth} \
            {params.paired} \
        2>> {log} 1>&2

        fa_files=$(find {output.bins} -name "*.fa")
        for fa in $fa_files; do
            pigz --best --verbose "$fa"
        done 2>> {log} 1>&2
        """


rule assemble__metabat2:
    """Run metabat2 over all assemblies"""
    input:
        [METABAT2 / assembly_id for assembly_id in ASSEMBLIES],
