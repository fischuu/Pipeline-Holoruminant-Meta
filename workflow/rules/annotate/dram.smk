rule _annotate__dram__annotate:
    """Annotate dereplicate genomes with DRAM"""
    input:
        dereplicated_genomes=DREP / "dereplicated_genomes",
        gtdbtk_summary=GTDBTK / "gtdbtk.summary.tsv",
        dram_db=features["databases"]["dram"],
    output:
        annotation=DRAM / "annotations.tsv",
        trnas=DRAM / "trnas.tsv",
        rrnas=DRAM / "rrnas.tsv",
        tarball=DRAM / "annotate.tar.gz",
    log:
        DRAM / "annotate.log",
    conda:
        "__environment__.yml"
    singularity:
        docker["annotate"]
    threads: 24
    params:
        min_contig_size=1500,
        out_dir=DRAM,
        tmp_dir=DRAM / "annotate",
        parallel_retries=5,
    resources:
        mem_mb=32 * 1024,
        runtime=48 * 60,
    shell:
        """
        rm \
            --recursive \
            --force \
            --verbose {params.tmp_dir} \
        2>> {log} 1>&2

        mkdir \
            --parents \
            {params.tmp_dir} \
        2>>{log} 1>&2

        DRAM-setup.py set_database_locations \
            --amg_database_loc          {input.dram_db}/amg_database.*.tsv \
            --dbcan_fam_activities_loc  {input.dram_db}/CAZyDB.*.fam-activities.txt \
            --dbcan_loc                 {input.dram_db}/dbCAN-HMMdb-V*.txt \
            --dbcan_subfam_ec_loc       {input.dram_db}/CAZyDB.*.fam.subfam.ec.txt \
            --description_db_loc        {input.dram_db}/description_db.sqlite \
            --etc_module_database_loc   {input.dram_db}/etc_mdoule_database.*.tsv \
            --function_heatmap_form_loc {input.dram_db}/function_heatmap_form.*.tsv \
            --genome_summary_form_loc   {input.dram_db}/genome_summary_form.*.tsv \
            --kofam_hmm_loc             {input.dram_db}/kofam_profiles.hmm \
            --kofam_ko_list_loc         {input.dram_db}/kofam_ko_list.tsv \
            --module_step_form_loc      {input.dram_db}/module_step_form.*.tsv \
            --peptidase_loc             {input.dram_db}/peptidases.*.mmsdb \
            --pfam_hmm_loc              {input.dram_db}/Pfam-A.hmm.dat.gz \
            --pfam_loc                  {input.dram_db}/pfam.mmspro \
            --viral_loc                 {input.dram_db}/refseq_viral.*.mmsdb \
            --vog_annotations_loc       {input.dram_db}/vog_annotations_latest.tsv.gz \
            --vogdb_loc                 {input.dram_db}/vog_latest_hmms.txt \
        2>> {log} 1>&2

        parallel \
            --jobs {threads} \
            --retries {params.parallel_retries} \
            DRAM.py annotate \
                --input_fasta {{}} \
                --output_dir {params.tmp_dir}/{{/.}} \
                --threads 1 \
                --gtdb_taxonomy {input.gtdbtk_summary} \
        ::: {input.dereplicated_genomes}/*.fa.gz \
        2>> {log} 1>&2

        for file in annotations trnas rrnas ; do
            ( csvstack \
                --tabs \
                {params.tmp_dir}/*/$file.tsv \
            | csvformat \
                --out-tabs \
            > {params.out_dir}/$file.tsv \
            ) 2>> {log}
        done

        tar \
            --create \
            --directory {params.out_dir} \
            --file {output.tarball} \
            --remove-files \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            annotate \
        2>> {log} 1>&2
        """


rule _annotate__dram__distill:
    """Distill DRAM annotations."""
    input:
        annotations=DRAM / "annotations.tsv",
        trnas=DRAM / "trnas.tsv",
        rrnas=DRAM / "rrnas.tsv",
        dram_db=features["databases"]["dram"],
    output:
        genome=DRAM / "genome_stats.tsv",
        metabolism=DRAM / "metabolism_summary.xlsx",
        product_html=DRAM / "product.html",
        product_tsv=DRAM / "product.tsv",
    log:
        DRAM / "distill.log2",
    conda:
        "__environment__.yml"
    singularity:
        docker["annotate"]
    resources:
        mem_mb=16 * 1024,
        runtime=24 * 60,
    params:
        outdir_tmp=DRAM / "distill",
        outdir=DRAM,
    shell:
        """
        DRAM-setup.py set_database_locations \
            --amg_database_loc          {input.dram_db}/amg_database.*.tsv \
            --dbcan_fam_activities_loc  {input.dram_db}/CAZyDB.*.fam-activities.txt \
            --dbcan_loc                 {input.dram_db}/dbCAN-HMMdb-V*.txt \
            --dbcan_subfam_ec_loc       {input.dram_db}/CAZyDB.*.fam.subfam.ec.txt \
            --description_db_loc        {input.dram_db}/description_db.sqlite \
            --etc_module_database_loc   {input.dram_db}/etc_mdoule_database.*.tsv \
            --function_heatmap_form_loc {input.dram_db}/function_heatmap_form.*.tsv \
            --genome_summary_form_loc   {input.dram_db}/genome_summary_form.*.tsv \
            --kofam_hmm_loc             {input.dram_db}/kofam_profiles.hmm \
            --kofam_ko_list_loc         {input.dram_db}/kofam_ko_list.tsv \
            --module_step_form_loc      {input.dram_db}/module_step_form.*.tsv \
            --peptidase_loc             {input.dram_db}/peptidases.*.mmsdb \
            --pfam_hmm_loc              {input.dram_db}/Pfam-A.hmm.dat.gz \
            --pfam_loc                  {input.dram_db}/pfam.mmspro \
            --viral_loc                 {input.dram_db}/refseq_viral.*.mmsdb \
            --vog_annotations_loc       {input.dram_db}/vog_annotations_latest.tsv.gz \
            --vogdb_loc                 {input.dram_db}/vog_latest_hmms.txt \
        2>> {log} 1>&2

        DRAM.py distill \
            --input_file {input.annotations} \
            --rrna_path {input.rrnas} \
            --trna_path {input.trnas} \
            --output_dir {params.outdir_tmp} \
        2> {log} 1>&2

        mv {params.outdir_tmp}/* {params.outdir}/ 2>> {log} 1>&2
        rmdir {params.outdir_tmp} 2>> {log} 1>&2
        """


rule annotate__dram:
    """Run DRAM on dereplicated genomes."""
    input:
        rules._annotate__dram__distill.output,
