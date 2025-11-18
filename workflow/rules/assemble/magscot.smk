rule assemble__magscot__prodigal:
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
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__magscot__prodigal")
    resources:
        runtime=esc("runtime", "assemble__magscot__prodigal"),
        mem_mb=esc("mem_mb", "assemble__magscot__prodigal"),
        cpus_per_task=esc("cpus", "assemble__magscot__prodigal"),
        slurm_partition=esc("partition", "assemble__magscot__prodigal"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__magscot__prodigal')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__magscot__prodigal"))
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


rule assemble__magscot__hmmsearch_pfam:
    """Run hmmsearch over the predicted proteins of an assembly using Pfam as database"""
    input:
        proteins=MAGSCOT / "{assembly_id}" / "prodigal.faa",
        hmm=features["magscot"]["pfam_hmm"],
    output:
        tblout=MAGSCOT / "{assembly_id}" / "pfam.tblout.gz",
    log:
        MAGSCOT / "{assembly_id}" / "pfam.log",
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__magscot__hmmsearch_pfam")
    resources:
        runtime=esc("runtime", "assemble__magscot__hmmsearch_pfam"),
        mem_mb=esc("mem_mb", "assemble__magscot__hmmsearch_pfam"),
        cpus_per_task=esc("cpus", "assemble__magscot__hmmsearch_pfam"),
        slurm_partition=esc("partition", "assemble__magscot__hmmsearch_pfam"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__magscot__hmmsearch_pfam')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__magscot__hmmsearch_pfam"))
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


rule assemble__magscot__hmmsearch_tigr:
    """Run hmmsearch over the predicted proteins of an assembly using TIGR as database"""
    input:
        proteins=MAGSCOT / "{assembly_id}" / "prodigal.faa",
        hmm=features["magscot"]["tigr_hmm"],
    output:
        tblout=MAGSCOT / "{assembly_id}" / "tigr.tblout.gz",
    log:
        MAGSCOT / "{assembly_id}" / "tigr.log",
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__magscot__hmmsearch_tigr")
    resources:
        runtime=esc("runtime", "assemble__magscot__hmmsearch_tigr"),
        mem_mb=esc("mem_mb", "assemble__magscot__hmmsearch_tigr"),
        cpus_per_task=esc("cpus", "assemble__magscot__hmmsearch_tigr"),
        slurm_partition=esc("partition", "assemble__magscot__hmmsearch_tigr"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__magscot__hmmsearch_tigr')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__magscot__hmmsearch_tigr"))
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


rule assemble__magscot__join_hmms:
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
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__magscot__join_hmms")
    resources:
        runtime=esc("runtime", "assemble__magscot__join_hmms"),
        mem_mb=esc("mem_mb", "assemble__magscot__join_hmms"),
        cpus_per_task=esc("cpus", "assemble__magscot__join_hmms"),
        slurm_partition=esc("partition", "assemble__magscot__join_hmms"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__magscot__join_hmms')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__magscot__join_hmms"))
    shell:
        """
        ( (zgrep -v "^#" {input.tigr_tblout} || true) \
        | awk '{{print $1 "\\t" $3 "\\t" $5}}' ) \
        >  {output.merged} 2>  {log}

        ( (zgrep -v "^#" {input.pfam_tblout} || true) \
        | awk '{{print $1 "\\t" $4 "\\t" $5}}' ) \
        >> {output.merged} 2>> {log}
        """


rule assemble__magscot__merge_contig_to_bin:
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
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__magscot__merge_contig_to_bin")
    resources:
        runtime=esc("runtime", "assemble__magscot__merge_contig_to_bin"),
        mem_mb=esc("mem_mb", "assemble__magscot__merge_contig_to_bin"),
        cpus_per_task=esc("cpus", "assemble__magscot__merge_contig_to_bin"),
        slurm_partition=esc("partition", "assemble__magscot__merge_contig_to_bin"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__magscot__merge_contig_to_bin')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__magscot__merge_contig_to_bin"))
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

rule assemble__magscot__run:
    """Run MAGSCOT over one assembly"""
    input:
        contigs_to_bin=MAGSCOT / "{assembly_id}" / "contigs_to_bin.tsv",
        hmm=MAGSCOT / "{assembly_id}" / "hmm.tblout",
    output:
        ar53=MAGSCOT / "{assembly_id}" / "magscot.gtdb_rel207_ar53.out",
        bac120=MAGSCOT / "{assembly_id}" / "magscot.gtdb_rel207_bac120.out",
        refined_contig_to_bin=MAGSCOT / "{assembly_id}" / "magscot.refined.contig_to_bin.out",
        refined_out=MAGSCOT / "{assembly_id}" / "magscot.refined.out",
        scores=MAGSCOT / "{assembly_id}" / "magscot.scores.out",
    log:
        MAGSCOT / "{assembly_id}/magscot.log",
    container:
        docker["assemble"]
    params:
        out_prefix=lambda w: MAGSCOT / w.assembly_id / "magscot",
        extra=params["assemble"]["magscot"]["extra"],
        th=params["assemble"]["magscot"]["threshold"],
        script_folder=SCRIPT_FOLDER,
    threads: esc("cpus", "assemble__magscot__run")
    resources:
        runtime=esc("runtime", "assemble__magscot__run"),
        mem_mb=esc("mem_mb", "assemble__magscot__run"),
        cpus_per_task=esc("cpus", "assemble__magscot__run"),
        slurm_partition=esc("partition", "assemble__magscot__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__magscot__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__magscot__run"))
    shell:
        """
        set -e
        {{
        Rscript --vanilla {params.script_folder}/MAGScoT/MAGScoT.R \
            --input {input.contigs_to_bin} \
            --hmm {input.hmm} \
            --out {params.out_prefix} \
            {params.extra} \
            --threshold {params.th} \
         2> {log} 1>&2 || echo "No result but proceeding."
         }}
        touch {output.ar53} {output.bac120} {output.refined_contig_to_bin} {output.refined_out} {output.scores};
        
        echo $? >> {log}
        
        """


rule assemble__magscot__reformat:
    """Reformat the results from MAGSCOT"""
    input:
        refined_contig_to_bin=MAGSCOT / "{assembly_id}" / "magscot.refined.contig_to_bin.out",
    output:
        clean=MAGSCOT / "{assembly_id}" / "magscot.reformat.tsv",
    log:
        MAGSCOT / "{assembly_id}" / "magscot.reformat.log",
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__magscot__reformat")
    resources:
        runtime=esc("runtime", "assemble__magscot__reformat"),
        mem_mb=esc("mem_mb", "assemble__magscot__reformat"),
        cpus_per_task=esc("cpus", "assemble__magscot__reformat"),
        slurm_partition=esc("partition", "assemble__magscot__reformat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__magscot__reformat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__magscot__reformat"))
    params:
        script_folder=SCRIPT_FOLDER,
    shell:
        """
        set -e
        
        if [ -s {input.refined_contig_to_bin} ]; then
            Rscript --vanilla {params.script_folder}/clean_magscot_bin_to_contig.R \
                --input-file {input.refined_contig_to_bin} \
               --output-file {output.clean} \
            2> {log} 1>&2
        else
            echo "Input file is empty or missing. Skipping reformat step." > {log}
            touch {output.clean}  # Create empty output file
        fi
        """


rule assemble__magscot__rename:
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
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__magscot__rename")
    resources:
        runtime=esc("runtime", "assemble__magscot__rename"),
        mem_mb=esc("mem_mb", "assemble__magscot__rename"),
        cpus_per_task=esc("cpus", "assemble__magscot__rename"),
        slurm_partition=esc("partition", "assemble__magscot__rename"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__magscot__rename')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__magscot__rename"))
    params:
        script_folder=SCRIPT_FOLDER,
    shell:
        """
        if [ -s {input.clean} ]; then
            python {params.script_folder}/reformat_fasta_magscot.py \
            <(gzip -dc {input.assembly}) \
            {input.clean} \
        | pigz \
            --best \
            > {output.fasta} 2>> {log}  # Ajoute les erreurs au log
        else
            echo "No data found, skipping renaming step." > {log}
            touch {output.fasta}  # Create an empty output file to avoid job failure
        fi
        """


rule assemble__magscot:
    """Run MAGSCOT over all assemblies"""
    input:
        [MAGSCOT / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
