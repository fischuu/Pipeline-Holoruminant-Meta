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
    ass_bowtie2 [label = 'Bowtie2']
    ass_concoct [label = 'Concoct']
    ass_drep [label = 'DRep']
    ass_magscot [label = 'MAGScot']
    ass_maxbin2 [label = 'MaxBin2']
    ass_megahit [label = 'MEGAHIT, v1.2.9']
    ass_metabat2 [label = 'MetaBAT2']
    ass_metaspades [label = 'MetaSPAdes, v4.0.0']

    # Define annotation nodes
    annotate_bakta [label = 'Bakta']
    annotate_quast [label = 'QUAST']
    annotate_gtdbtk [label = 'GTDB-Tk']
    annotate_dram [label = 'DRAM']
    annotate_eggnog [label = 'EggNOG']
    annotate_checkm2 [label = 'CheckM2']
    annotate_proteinortho [label = 'Proteinortho']
    annotate_phylophlan [label = 'PhyloPhlan']

    # Define a representative node for the user input boxes
    pre_host [label = 'Host genomes', style=dotted, color=gray]
    ass_assembler_rep [label = 'Select Assembler', style=dotted, color=gray]

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
    
    # Define edges for assembly
    pre_bowtie2 -> ass_assembler_rep
    ass_bowtie2 -> ass_concoct
    ass_concoct -> ass_drep
    ass_drep -> ass_magscot
    ass_magscot -> ass_maxbin2
    ass_maxbin2 -> ass_megahit
    ass_megahit -> ass_metabat2
    ass_metabat2 -> ass_metaspades

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

      # Define a subgraph to group assembler nodes
      subgraph cluster_assembler {
        label = 'Assembler'
        style = dashed
        color = lightgray

        ass_megahit
        ass_metaspades
      }
    }

    # Define a subgraph to group reads nodes
    subgraph cluster_reads {
      label = 'Reads'
      style = dashed
      color = lightgray

      link
      fastqc
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

    # Define edge from DRep to Annotate
    ass_drep -> cluster_annotate
  }
")

# Render the graph
graph
