rule assemble__drep__separate_bins:
    input:
        assemblies=[MAGSCOT / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
    output:
        out_dir=directory(DREP / "separated_bins"),
    log:
        DREP / "separate_bins.log",
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__drep__separate_bins")
    resources:
        runtime=esc("runtime", "assemble__drep__separate_bins"),
        mem_mb=esc("mem_mb", "assemble__drep__separate_bins"),
        cpus_per_task=esc("cpus", "assemble__drep__separate_bins"),
        slurm_partition=esc("partition", "assemble__drep__separate_bins"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'assemble__drep__separate_bins')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__drep__separate_bins"))
    shell:
        """
        mkdir --parents {output.out_dir} 2> {log} 1>&2

        ( gzip \
            --decompress \
            --stdout \
            {input.assemblies} \
        | paste - - \
        | tr -d ">" \
        | tr "@" "\t" \
        | awk \
            '{{print ">" $1 "@" $2 "\\n" $3 > "{output.out_dir}/" $1 ".fa" }}' \
        ) >> {log} 2>&1
        """


rule assemble__drep__run:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=DREP / "separated_bins",
    output:
        dereplicated_genomes=directory(DREP / "dereplicated_genomes"),
        data=DREP / "data.tar.gz",
        data_tables=DREP / "data_tables.tar.gz",
    log:
        DREP / "drep.log",
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__drep__run")
    resources:
        runtime=esc("runtime", "assemble__drep__run"),
        mem_mb=esc("mem_mb", "assemble__drep__run"),
        cpus_per_task=esc("cpus", "assemble__drep__run"),
        slurm_partition=esc("partition", "assemble__drep__run"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'assemble__drep__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__drep__run"))
    params:
        out_dir=DREP,
        completeness=params["assemble"]["drep"]["completeness"],
        contamination=params["assemble"]["drep"]["contamination"],
        P_ani=params["assemble"]["drep"]["P_ani"],
        S_ani=params["assemble"]["drep"]["S_ani"],
        nc=params["assemble"]["drep"]["nc"],
        extra=params["assemble"]["drep"]["extra"]
    shell:
        """
        
        rm \
            --recursive \
            --force \
            {params.out_dir}/data_tables \
            {params.out_dir}/data \
            {params.out_dir}/dereplicated_genomes \
            {params.out_dir}/figures \
            {params.out_dir}/log \
        2> {log}.{resources.attempt} 1>&2

        echo "This is the current TMPDIR: " 2>> {log}.{resources.attempt} 1>&2
        echo $TMPDIR 2>> {log}.{resources.attempt} 1>&2
        
        echo "Value of LOCAL_SCRATCH before export: $LOCAL_SCRATCH" 2>> {log}.{resources.attempt} 1>&2

        export TMPDIR=$LOCAL_SCRATCH 2>> {log}.{resources.attempt} 1>&2
        
        echo "This is the current TMPDIR: " 2>> {log}.{resources.attempt} 1>&2
        echo $TMPDIR 2>> {log}.{resources.attempt} 1>&2
        
        dRep dereplicate \
            {params.out_dir} \
            --processors {threads} \
            --completeness {params.completeness} \
            --contamination {params.contamination} \
            --P_ani {params.P_ani} \
            --S_ani {params.S_ani} \
            -nc {params.nc} \
            {params.extra} \
            --genomes {input.genomes}/*.fa \
        2>> {log}.{resources.attempt} 1>&2

        for folder in data data_tables ; do
            tar \
                --create \
                --directory {params.out_dir} \
                --file {params.out_dir}/${{folder}}.tar.gz \
                --remove-files \
                --use-compress-program="pigz --processes {threads}" \
                --verbose \
                ${{folder}} \
            2>> {log} 1>&2
        done

        files=$(find {output.dereplicated_genomes} -type f ! -name "*.gz")
        for file in $files; do
            gzip "$file"
        done

        mv {log}.{resources.attempt} {log}
        """


rule assemble__drep__join_genomes:
    """Join all the dereplicated genomes into a single file."""
    input:
        DREP / "dereplicated_genomes",
    output:
        DREP / "dereplicated_genomes.fa.gz",
    log:
        DREP / "dereplicated_genomes.log",
    container:
        docker["assemble"]
    threads: esc("cpus", "assemble__drep__join_genomes")
    resources:
        runtime=esc("runtime", "assemble__drep__join_genomes"),
        mem_mb=esc("mem_mb", "assemble__drep__join_genomes"),
        cpus_per_task=esc("cpus", "assemble__drep__join_genomes"),
        slurm_partition=esc("partition", "assemble__drep__join_genomes"),
        slurm_extra=lambda wc, attempt: f"--gres=nvme:{get_resources(wc, attempt, 'assemble__drep__join_genomes')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__drep__join_genomes"))
    shell:
        """
        ( zcat \
            {input}/*.fa.gz \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output} \
        ) 2> {log}
        """


rule assemble__drep:
    input:
        DREP / "dereplicated_genomes.fa.gz",
