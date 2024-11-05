rule _contig_annotate__cramToBam_:
    """Create temporary bam-files for quantification"""
    input:
        ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram"
    output:
        temp(ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam"),
    log:
        log=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cramToBam.log",
    conda:
        "__environment__.yml"
    container:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["single_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["single_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        samtools index {input} 2> {log} 1>&2
        samtools view -b -o {output} {input} 2> {log} 1>&2
        """


rule _contigAnnotate__featureCounts:
    """
    Quantify the mapped assembly reads in predicted genes (featureCounts).
    """
    input:
        bam=ASSEMBLE_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam",
        gtf=CONTIG_PRODIGAL / "{assembly_id}/{assembly_id}.prodigal.gtf",
    output:
        file=CONTIG_FEATURECOUNTS / "{assembly_id}/{assembly_id}_{sample_id}.{library_id}.prodigal_fc.txt",
    log:
        CONTIG_FEATURECOUNTS / "{assembly_id}/{assembly_id}_{sample_id}.{library_id}.prodigal_fc.log",
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
    container: docker["quantify"]
    shell:"""
        featureCounts -p \
                      -T {threads} \
                      -a {input.gtf} \
                      -o {output.file} \
                      -t CDS \
                      -g ID \
                      {input.bam} 2> {log}
    """
    
    
rule contig_annotate__featurecount:
    """Quantify all samples for all the assemblies that they belong to"""
    input:
        [
            CONTIG_FEATURECOUNTS / "{assembly_id}/{assembly_id}_{sample_id}.{library_id}.prodigal_fc.txt"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
