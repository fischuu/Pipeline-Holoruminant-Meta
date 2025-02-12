rule _assemble__magscot__prodigal:
    """Run prodigal over a single assembly"""
    input:
        assembly=lambda wildcards: (
            MEGAHIT / f"{wildcards.assembly_id}.fa.gz" if config["assembler"] == "megahit" else
            METASPADES / f"{wildcards.assembly_id}.fa.gz"
        ),
    output:
        proteins=MAGSCOT / "{assembly_id}" / "prodigal.faa",
    log:
        MAGSCOT / "{assembly_id}" / "prodigal.log",
    conda:
        "__environment__.yml"
    container:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        nvme = config["resources"]["nvme"]["small"],
        attempt=get_attempt,
    retries: 5
    shell:
        """
        
        ( gzip \
            --decompress \
            --stdout \
            {input.assembly} \
        | parallel \
            --jobs {threads} \
            --block 1M \
            --recstart '>' \
            --pipe \
            --keep-order \
            prodigal \
                -p meta \
                -a /dev/stdout \
                -d /dev/null  \
                -o /dev/null \
        > {output.proteins} \
        ) 2> {log}.{resources.attempt}

        mv {log}.{resources.attempt} {log}
        """


rule _assemble__magscot__hmmsearch_pfam:
    """Run hmmsearch over the predicted proteins of an assembly using Pfam as database"""
    input:
        proteins=MAGSCOT / "{assembly_id}" / "prodigal.faa",
        hmm=features["magscot"]["pfam_hmm"],
    output:
        tblout=MAGSCOT / "{assembly_id}" / "pfam.tblout.gz",
    log:
        MAGSCOT / "{assembly_id}" / "pfam.log",
    conda:
        "__environment__.yml"
    container:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        hmmsearch \
            -o /dev/null \
            --tblout >(pigz --best --processes {threads} > {output.tblout}) \
            --noali \
            --notextw \
            --cut_nc \
            --cpu {threads} \
            {input.hmm} \
            {input.proteins} \
        2> {log} 1>&2
        """


rule _assemble__magscot__hmmsearch_tigr:
    """Run hmmsearch over the predicted proteins of an assembly using TIGR as database"""
    input:
        proteins=MAGSCOT / "{assembly_id}" / "prodigal.faa",
        hmm=features["magscot"]["tigr_hmm"],
    output:
        tblout=MAGSCOT / "{assembly_id}" / "tigr.tblout.gz",
    log:
        MAGSCOT / "{assembly_id}" / "tigr.log",
    conda:
        "__environment__.yml"
    container:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        hmmsearch \
            -o /dev/null \
            --tblout >(pigz --best --processes {threads} > {output.tblout}) \
            --noali \
            --notextw \
            --cut_nc \
            --cpu {threads} \
            {input.hmm} \
            {input.proteins} \
        2> {log} 1>&2
        """


rule _assemble__magscot__join_hmms:
    """Join the results of hmmsearch over TIGR and Pfam

    Note: "|| true" is used to avoid grep returning an error code when no lines are found
    """
    input:
        tigr_tblout=MAGSCOT / "{assembly_id}" / "tigr.tblout.gz",
        pfam_tblout=MAGSCOT / "{assembly_id}" / "pfam.tblout.gz",
    output:
        merged=MAGSCOT / "{assembly_id}" / "hmm.tblout",
    log:
        MAGSCOT / "{assembly_id}" / "hmm.log",
    conda:
        "__environment__.yml"
    container:
        docker["assemble"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        ( (zgrep -v "^#" {input.tigr_tblout} || true) \
        | awk '{{print $1 "\\t" $3 "\\t" $5}}' ) \
        >  {output.merged} 2>  {log}

        ( (zgrep -v "^#" {input.pfam_tblout} || true) \
        | awk '{{print $1 "\\t" $4 "\\t" $5}}' ) \
        >> {output.merged} 2>> {log}
        """


rule _assemble__magscot__merge_contig_to_bin:
    """Merge the contig to bin files from CONCOCT, MaxBin2 and MetaBAT2

    The output file should have the following format:
    BIN_ID <TAB> CONTIG_ID <TAB> METHOD
    """
    input:
        concoct=CONCOCT / "{assembly_id}",
        maxbin2=MAXBIN2 / "{assembly_id}",
        metabat2=METABAT2 / "{assembly_id}",
    output:
        MAGSCOT / "{assembly_id}" / "contigs_to_bin.tsv",
    log:
        MAGSCOT / "{assembly_id}" / "contigs_to_bin.log",
    conda:
        "__environment__.yml"
    container:
        docker["assemble"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        for file in $(find {input.concoct} -name "*.fa.gz" -type f) ; do
            bin_id=$(basename $file .fa)
            zgrep ^">" $file | tr -d ">" \
            | awk -v bin_id=$bin_id '{{print "bin_" bin_id "\\t" $1 "\\tconcoct"}}'
        done > {output} 2> {log}

        for file in $(find {input.maxbin2} -name "*.fa.gz" -type f) ; do
            bin_id=$(basename $file .fa)
            zgrep ^">" $file | tr -d ">" \
            | awk -v bin_id=$bin_id '{{print "bin_" bin_id "\\t" $1 "\\tmaxbin2"}}'
        done >> {output} 2>> {log}

        for file in $(find {input.metabat2} -name "*.fa.gz" -type f) ; do
            bin_id=$(basename $file .fa)
            zgrep ^">" $file | tr -d ">" \
            | awk -v bin_id=$bin_id '{{print "bin_" bin_id "\\t" $1 "\\tmetabat2"}}'
        done >> {output} 2>> {log}
        """

rule _assemble__magscot__run:
    """Run MAGSCOT over one assembly"""
    input:
        contigs_to_bin=MAGSCOT / "{assembly_id}" / "contigs_to_bin.tsv",
        hmm=MAGSCOT / "{assembly_id}" / "hmm.tblout",
    output:
        ar53=MAGSCOT / "{assembly_id}" / "magscot.gtdb_rel207_ar53.out",
        bac120=MAGSCOT / "{assembly_id}" / "magscot.gtdb_rel207_bac120.out",
        refined_contig_to_bin=MAGSCOT
        / "{assembly_id}"
        / "magscot.refined.contig_to_bin.out",
        refined_out=MAGSCOT / "{assembly_id}" / "magscot.refined.out",
        scores=MAGSCOT / "{assembly_id}" / "magscot.scores.out",
    log:
        MAGSCOT / "{assembly_id}/magscot.log",
    conda:
        "__environment__.yml"
    container:
        docker["assemble"]
    params:
        out_prefix=lambda w: MAGSCOT / w.assembly_id / "magscot",
        extra=params["assemble"]["magscot"]["extra"],
        th=params["assemble"]["magscot"]["threshold"],
        script_folder=SCRIPT_FOLDER,
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"],
        time =  config["resources"]["time"]["longrun"],
    shell:
        """
        set -e
        
        Rscript --vanilla {params.script_folder}/MAGScoT/MAGScoT.R \
            --input {input.contigs_to_bin} \
            --hmm {input.hmm} \
            --out {params.out_prefix} \
            {params.extra} \
            --threshold {params.th} \
         2> {log} 1>&2 
         
        echo $? >> {log}
        
        """


rule _assemble__magscot__reformat:
    """Reformat the results from MAGSCOT"""
    input:
        refined_contig_to_bin=MAGSCOT
        / "{assembly_id}"
        / "magscot.refined.contig_to_bin.out",
    output:
        clean=MAGSCOT / "{assembly_id}" / "magscot.reformat.tsv",
    log:
        MAGSCOT / "{assembly_id}" / "magscot.reformat.log",
    conda:
        "__environment__.yml"
    container:
        docker["assemble"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["lowmem"],
        time =  config["resources"]["time"]["shortrun"],
    params:
        script_folder=SCRIPT_FOLDER,
    shell:
        """
        set -e
        
        Rscript --vanilla {params.script_folder}/clean_magscot_bin_to_contig.R \
            --input-file {input.refined_contig_to_bin} \
            --output-file {output.clean} \
        2> {log} 1>&2
        """


rule _assemble__magscot__rename:
    """Rename the contigs in the assembly to match the assembly and bin names"""
    input:
        assembly=lambda wildcards: (
            MEGAHIT / f"{wildcards.assembly_id}.fa.gz" if config["assembler"] == "megahit" else
            METASPADES / f"{wildcards.assembly_id}.fa.gz"
        ),
        clean=MAGSCOT / "{assembly_id}" / "magscot.reformat.tsv",
    output:
        fasta=MAGSCOT / "{assembly_id}.fa.gz",
    log:
        MAGSCOT / "{assembly_id}" / "magscot.rename.log",
    conda:
        "__environment__.yml"
    container:
        docker["assemble"]
    resources:
        mem_per_cpu=config["resources"]["mem_per_cpu"]["lowmem"],
        time =  config["resources"]["time"]["shortrun"],
    params:
        script_folder=SCRIPT_FOLDER,
    shell:
        """
        ( python {params.script_folder}/reformat_fasta_magscot.py \
            <(gzip -dc {input.assembly}) \
            {input.clean} \
        | pigz \
            --best \
        > {output.fasta} \
        ) 2> {log}
        """


rule assemble__magscot:
    """Run MAGSCOT over all assemblies"""
    input:
        [MAGSCOT / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
