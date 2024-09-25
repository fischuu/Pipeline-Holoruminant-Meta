rule _contig_annotate__eggnog_find_homology:
    """
    Find homolog genes in data (EGGNOG).
    """
    input:
        folder = CONTIG_PRODIGAL / "{assembly_id}/Chunks/",
        files = CONTIG_PRODIGAL / "{assembly_id}/Chunks/prodigal.chunk.{i}"
    output:
        file=CONTIG_EGGNOG / "{assembly_id}/Chunks/prodigal.chunk.{i}.emapper.seed_orthologs"
    log:
        CONTIG_EGGNOG / "{assembly_id}/Chunks/prodigal.chunk.{i}.log"
    params:
        folder = lambda wildcards: CONTIG_EGGNOG / f"{wildcards.assembly_id}/Chunks/",
        tmp=config["nvme_storage"],
        out="prodigal.chunk.{i}",
        fa=features["databases"]["eggnog"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
        nvme = config["resources"]["nvme"]["large"]
    container:
        docker["annotate"]
    shell:""" 
         cp {params.fa}/eggnog* {params.tmp}  &> {log};
         emapper.py -m diamond --data_dir {params.tmp} --no_annot --no_file_comments --cpu {threads} -i {input.files} --output_dir {params.folder} -o {params.out}  2>> {log} 1>&2;
    """
    
                  
rule _contig_annotate__aggregate_assemblies_eggnog:
    input:
        _contig_annotate_aggregate_assembly_eggnog_search,
    output:
        CONTIG_EGGNOG / "{assembly_id}/prodigal.emapper.seed_orthologs"
    log:
        CONTIG_EGGNOG / "{assembly_id}/prodigal.emapper.seed_orthologs.log"
    container:
        docker["annotate"]
    shell:"""
       cat {input} > {output} 2> {log}
    """
    
rule _contig_annotate__eggnog_orthology:
    """
    Find orthology and annotate (EGGNOG).
    """
    input:
        CONTIG_EGGNOG / "{assembly_id}/prodigal.emapper.seed_orthologs",
    output:
        CONTIG_EGGNOG / "{assembly_id}/eggnog_output.emapper.annotations",
    log:
        CONTIG_EGGNOG / "{assembly_id}/prodigal.emapper.annotations.log",
    params:
        tmp=config["nvme_storage"],
        fa=features["databases"]["eggnog"],
        out = lambda wildcards: CONTIG_EGGNOG / f"{wildcards.assembly_id}/eggnog_output",
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
        nvme = config["resources"]["nvme"]["large"]
    container:
        docker["annotate"]
    shell:"""
    # Actually, not sure if that makes any sense, I think the files should not be concatenated here, but treated still here concatenated and then rather be merged afterwards
    # For now it is alright, as we annotated approx 650 per seconds (=run takes around a day), but if you ever rerun this step again, change it that it will be processed on the chunks
    # and then merge the chunks!
       cp {params.fa}/eggnog* {params.tmp}  &>> {log};
       emapper.py --data_dir {params.tmp} --annotate_hits_table {input} --no_file_comments -o {params.out} --cpu {threads} &> {log}
    """
    ß
rule contig_assemble__eggnog:
    """Run eggnog on all assemblies"""
    input:
        [CONTIG_EGGNOG / "{assembly_id}/eggnog_output.emapper.annotations" for assembly_id in ASSEMBLIES],
