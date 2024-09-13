library(DiagrammeR)

# Create a DiagrammeR graph
graph <- grViz("
  digraph snakemake_pipeline {
    graph [layout = dot, rankdir = LR]

    node [shape = box, style = filled, fillcolor = lightblue]
    edge [arrowhead = vee]

    # Define new reads nodes
    input_reads [label = 'Input Reads', shape = folder]
    link [label = 'Link \n Assign names']
    fastqc [label = 'FastQC, v0.12.1 \n Raw reads QC']

    # Define preprocessing nodes
    pre_bowtie2 [label = 'Bowtie2, v2.5.1 \n Decontamination']
    pre_fastp [label = 'Fastp, v0.23.4 \n QC trimming']
    pre_fastqc [label = 'FastQC, v0.12.1']
    pre_kraken2 [label = 'Kraken2, v2.1.3 \n RefSeqV205_Complete \n Standard_20240112']
    pre_humann [label = 'HumanN3, v3.9 \n Chocophlan v201901_v31 \n Uniref uniref90_201901b_full']
    pre_metaphlan [label = 'MetaPhlAn, v4.1.1 \n mpa_vJun23_CHOCOPhlAnSGB_202307']
    pre_phyloflash [label = 'PhyloFlash, v3.4.2 \n SILVA_SSU.noLSU 138.1']
    pre_nonpareil [label = 'Nonpareil, v.3.4.1']
    pre_singlem [label = 'SingleM, v0.18.0 \n S4.3.0.GTDB_r220.metapackage_20240523']

    # Define assembly nodes
    ass_bowtie2 [label = 'Bowtie2, v2.5.1']
    ass_concoct [label = 'Concoct, v1.1.0ß']
    ass_drep [label = 'DRep, v3.4.3']
    ass_magscot [label = 'MAGScot, v???']
    ass_maxbin2 [label = 'MaxBin2, v2.2.7']
    ass_megahit [label = 'MEGAHIT, v1.2.9']
    ass_metabat2 [label = 'MetaBAT2, v2.15']
    ass_metaspades [label = 'MetaSPAdes, v4.0.0']

    # Define annotation nodes
    annotate_bakta [label = 'Bakta, v1.9.3', fillcolor = 'red']
    annotate_quast [label = 'QUAST, v5.2.0']
    annotate_gtdbtk [label = 'GTDB-Tk, v2.4.0 \n Release 220']
    annotate_dram [label = 'DRAM, v1.5.0 \n 20240524', fillcolor = 'red']
    annotate_eggnog [label = 'EggNOG, v2.1.12']
    annotate_checkm2 [label = 'CheckM2, v1.0.2 \n uniref100.KO.1.dmnd']
    annotate_proteinortho [label = 'Proteinortho, v6.3.1', fillcolor = 'red']
    annotate_phylophlan [label = 'PhyloPhlan, v3.1.1', fillcolor = 'red']

    # Define the contig annotate nodes
    contig_annotate_prodigal [label = 'Prodigal, v2.6.3']

#    dram_dbs [label = 'DRAM DBs \n
#    kegg: null \n
#    kofam_hmm: 20240524/kofam_profiles.hmm \n
#    kofam_ko_list: 20240524/kofam_ko_list.tsv \n
#    uniref: 20240524/uniref90.20240524.mmsdb \n
#    pfam: 20240524/pfam.mmspro \n
#    dbcan: 20240524/dbCAN-HMMdb-V11.txt \n
#    viral: 20240524/refseq_viral.20240524.mmsdb \n
#    peptidase: 20240524/peptidases.20240524.mmsdb \n
#    vogdb: 20240524/vog_latest_hmms.txt \n
#    camper_hmm: null \n
#    camper_fa_db: null \n
#    camper_hmm_cutoffs: null \n
#    camper_fa_db_cutoffs: null
#    ', shape = folder]
    
    # Define quantify nodes
    quantify_bowtie2 [label = 'Bowtie2, v2.5.1']
    quantify_coverm [label = 'CoverM, v0.6.1']
    quantify_samtools [label = 'Samtools, v1.18']

    # Define feature request nodes
    feature_phylophlan [label = 'Phylophlan', fillcolor = 'yellow']
    feature_dram2 [label = 'DRAM2', fillcolor = 'yellow']

    # Define a representative node for the user input boxes
    pre_host [label = 'Host genomes', style=dotted, color=gray]
    ass_assembler_rep [label = 'Select Assembler', style=dotted, color=gray]
    ass_assembly [label = 'Metagenome Assembly', style=dotted, color=gray]


    # Define edges for input reads
    input_reads -> link
    input_reads -> fastqc
    link -> pre_fastp
    

    # Define edges for preprocessing
    pre_fastp -> pre_host
    pre_host -> pre_bowtie2
    pre_fastp -> pre_fastqc
    pre_bowtie2 -> pre_fastqc
    pre_fastp -> pre_kraken2
    pre_bowtie2 -> pre_humann
    pre_bowtie2 -> pre_metaphlan
    pre_bowtie2 -> pre_phyloflash
    pre_bowtie2 -> pre_nonpareil
    pre_bowtie2 -> pre_singlem
    pre_bowtie2 -> ass_assembler_rep
    
    # Define edges for assembly
    ass_assembly -> ass_bowtie2
    ass_assembly -> ass_concoct
    ass_bowtie2 -> ass_concoct
    ass_concoct -> ass_magscot
    ass_metabat2 -> ass_magscot
    ass_maxbin2 -> ass_magscot
    ass_magscot -> ass_drep
    ass_assembly -> ass_magscot
    ass_assembly -> ass_maxbin2
    ass_bowtie2 -> ass_maxbin2
    ass_megahit -> ass_assembly
    ass_metaspades -> ass_assembly
    ass_assembly -> ass_metabat2
    ass_bowtie2 -> ass_metabat2
    pre_bowtie2 -> ass_bowtie2

    # Define edges for quantify
    ass_drep -> quantify_bowtie2
    ass_drep -> quantify_coverm
    quantify_bowtie2 -> quantify_coverm
    ass_drep -> quantify_samtools
    quantify_bowtie2 -> quantify_samtools
    
    # Define the edges for annotate
    ass_drep -> annotate_bakta
    ass_drep -> annotate_checkm2
    annotate_gtdbtk -> annotate_dram
    ass_drep -> annotate_dram
    ass_drep -> annotate_eggnog
    ass_drep -> annotate_gtdbtk
    ass_drep -> annotate_phylophlan
    annotate_bakta -> annotate_proteinortho
    ass_drep -> annotate_quast
#    dram_dbs -> annotate_dram

    # Define the edges for contig annotate
    ass_assembly -> contig_annotate_prodigal

    # Define edges between nodes and the representative node
    ass_assembler_rep -> ass_megahit
    ass_assembler_rep -> ass_metaspades

    # Define a subgraph to group preprocessing nodes
    subgraph cluster_preprocess {
      label = 'Preprocess, v0.5'
      style = dashed
      color = lightgray

      pre_bowtie2
      pre_fastp
      pre_fastqc
      pre_kraken2
      pre_humann
      pre_metaphlan
      pre_phyloflash
      pre_nonpareil
      pre_singlem
      pre_host
    }

    # Define a subgraph to group assembly nodes
    subgraph cluster_assemble {
      label = 'Assemble, v0.3'
      style = dashed
      color = lightgray

      ass_bowtie2
      ass_concoct
      ass_drep
      ass_magscot
      ass_maxbin2
      ass_assembly
      ass_metabat2

      # Define a subgraph to group assembler nodes
      subgraph cluster_assembler {
        label = 'Assembler'
        style = dashed
        color = lightgray

        ass_megahit
        ass_metaspades
      }
      
       # Define a subgraph to group binner nodes
      subgraph cluster_binner {
        label = 'Binner'
        style = dashed
        color = lightgray

        ass_concoct
        ass_metabat2
        ass_maxbin2
      }
    }

# Define a subgraph to group quantify nodes
    subgraph cluster_quantify {
      label = 'Quantify, v'
      style = dashed
      color = lightgray

      quantify_bowtie2
      quantify_coverm
      quantify_samtools
    }

    # Define a subgraph to group reads nodes
    subgraph cluster_reads {
      label = 'Reads'
      style = dashed
      color = lightgray

      link
      fastqc
    }

    # Define a subgraph to group contig annotate nodes
    subgraph cluster_contig_annotate {
      label = 'Contig annotate'
      style = dashed
      color = lightgray

      contig_annotate_prodigal
    }


    # Define a subgraph to group feature request nodes
    subgraph cluster_featureRequest {
      label = 'Feature Requests'
      style = dashed
      color = lightgray

      feature_dram2
      feature_phylophlan
    }

    # Define a subgraph to group annotation nodes
    subgraph cluster_annotate {
      label = 'Annotate'
      style = dashed
      color = lightgray

      annotate_bakta
      annotate_quast
      annotate_gtdbtk
      annotate_dram
      annotate_eggnog
      annotate_checkm2
      annotate_proteinortho
      annotate_phylophlan
    }
  }
")

# Render the graph
graph
