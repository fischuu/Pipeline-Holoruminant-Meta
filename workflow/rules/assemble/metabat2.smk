rule _assemble__metabat2__run:
    """Run metabat2 end-to-end on a single assembly"""
    input:
        crams=get_crams_from_assembly_id,
        crais=get_crais_from_assembly_id,
        assembly=MEGAHIT / "{assembly_id}.fa.gz",
    output:
        bins=directory(METABAT2 / "{assembly_id}"),
    log:
        METABAT2 / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["assemble"]
    params:
        bins_prefix=lambda w: METABAT2 / f"{w.assembly_id}/bin",
        bams=compose_bams_for_metabat2_run,
        depth=lambda w: METABAT2 / f"{w.assembly_id}.depth",
        paired=lambda w: METABAT2 / f"{w.assembly_id}.paired",
        workdir=METABAT2,
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
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
            --verbose \
        2> {log} 1>&2

        rm \
            --force \
            --verbose \
            {params.bams} \
            {params.depth} \
            {params.paired} \
        2>> {log} 1>&2

        find \
            {output.bins} \
            -name "*.fa" \
            -exec pigz --best --verbose {{}} \
        2>> {log} 1>&2
        """


rule assemble__metabat2:
    """Run metabat2 over all assemblies"""
    input:
        [METABAT2 / assembly_id for assembly_id in ASSEMBLIES],
