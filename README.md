# Overview
This Snakemake pipeline dedicated to Metagenomic data analysis consists out of several modules that cover a) read-based b) contig-based and c) MAG-based analyses as well as quantification, quality checks and a reporting module. Naturally, it runs seamlessly on HPC systems and all required software tools are bundled in docker container and/or conda environments. Further, all required databases are pre-configured and ready to be downloaded from a central place.

![alt text](https://github.com/fischuu/Pipeline-Holoruminant-Meta/blob/main/flowchart/flowchart.png?raw=true)

(Red marked rules have currently still unsolved issues)

# Requirements
The pipeline requires version 8 or later (Snakemake > 8.x)

Supports:
SLURM executor / local execution
conda environment (not tested)
docker/singularity/apptainer support

## Python dependencies
Since Snakemake 8, it is required to install a cluster-generic plugin to submit jobs to a queueing system of a HPC system. Please ensure you have installed the corresponding Snakemake plugin installed in case you want to submit your jobs to a queueing system

```
pip install snakemake-executor-plugin-cluster-generic
```

Of course you can also install specific plugins like the slurm plugin, but this might need more adjustments to the existing files.

Depending on your Python version, you need to install a Pandas version > 2.1, there were errors when newer Python versions met older Pandas version. In case you run into obscure Pandas error, please make sure to install a newer pandas, e.g.

```
pip install pandas==2.2.3
```

# Installation

You can install the pipeline by cloning this repository

The recommended setup is to have a dedicated pipeline folder (the cloned repository), that carries the functionality and which should not require any changes. 

Then the project should have somewhere an own folder and the required configuration files are copied to it. The steps to perform are

```
# Go to the folder, to where you would like to clone the pipeline, e.g. 
 cd /users/fischerd/git

# First, clone the pipeline into that folder
  git clone git@github.com:fischuu/Pipeline-Holoruminant-Meta.git

# In case the previous steps fails with an error that contains
# git@github.com: Permission denied (publickey)
# it indicates that you do not have a ssh key exchanged with GitHub and you could clone the repository then instead like this
# git clone https://github.com/fischuu/Pipeline-Holoruminant-Meta.git
  
# Setting ENV variable to get downstream code more generic (so, this is the directory to where you cloned the pipeline)
  cd Pipeline-Holoruminant-Meta
  PIPELINEFOLDER=$(pwd)
  
# If previous doesn't work, you can set it also manually like for example this
  PIPELINEFOLDER="/users/fischerd/git/Pipeline-Holoruminant-Meta"
```

Next, we setup a project folder in our scratch space of the HPC, here we will run the pipeline

```
# Go to the project space of your HPC, e.g. 
  cd /scratch/project_2009831

# Create a folder for the new project
  mkdir My_holor_project
  cd My_holor_project   
  
# For convenience, we set again a ENV variable, so that the code will be more generic
  PROJECTFOLDER=$(pwd)
  
# Or manually the same thing:  
  PROJECTFOLDER="/scratch/project_2009831/My_holor_project"
```

Then we need to download the precompiled databases and reference genomes.
Be prepared that this step will take some time (3 days) and disc space (3TB).
In case you have quick, local nvme discs, it is advisable to use them for
unpacking the files, as this will significantly increase the speed.

```
# Change to the project folder and prepare folders
  cd $PROJECTFOLDER
  mkdir -p reads
  mkdir -p resources/databases
  mkdir -p resources/reference

# Download the various pre-prepared reference databases
  cd $PROJECTFOLDER/resources/databases
  wget https://a3s.fi/Holoruminant-data/2024.09.18.bakta.tar.gz
  wget https://a3s.fi/Holoruminant-data/2024.09.18.diamond.tar.gz
  wget https://a3s.fi/Holoruminant-data/2024.09.18.eggnog.tar.gz
  wget https://a3s.fi/Holoruminant-data/2024.09.18.humann.tar.gz
  wget https://a3s.fi/Holoruminant-data/2024.09.18.metaphlan4.tar.gz
  wget https://a3s.fi/Holoruminant-data/2024.09.18.phylophlan.tar.gz
  wget https://a3s.fi/Holoruminant-data/2024.09.18.checkm2.tar.gz
  wget https://a3s.fi/Holoruminant-data/2025.01.28.dram.tar.gz
  wget https://a3s.fi/Holoruminant-data/2024.09.18.gtdbtk.tar.gz
  wget https://a3s.fi/Holoruminant-data/2024.09.18.kraken2.tar.gz
  wget https://a3s.fi/Holoruminant-data/2025.02.14.krona.tar.gz
  wget https://a3s.fi/Holoruminant-data/2024.09.18.phyloflash.tar.gz
  wget https://a3s.fi/Holoruminant-data/2024.09.18.singlem.tar.gz

# Unpack all the databases
  tar -xvf 2024.09.18.bakta.tar.gz
  tar -xvf 2024.09.18.diamond.tar.gz
  tar -xvf 2024.09.18.eggnog.tar.gz
  tar -xvf 2024.09.18.humann.tar.gz
  tar -xvf 2024.09.18.metaphlan4.tar.gz
  tar -xvf 2024.09.18.phylophlan.tar.gz
  tar -xvf 2024.09.18.checkm2.tar.gz
  tar -xvf 2025.01.28.dram.tar.gz
  tar -xvf 2024.09.18.gtdbtk.tar.gz
  tar -xvf 2024.09.18.kraken2.tar.gz
  tar -xvf 2025.02.14.krona.tar.gz
  tar -xvf 2024.09.18.phyloflash.tar.gz
  tar -xvf 2024.09.18.singlem.tar.gz

# Get the reference genomes relevant for Holorumiant for host contamination removal
# Obviously, you can also use your own set of reference genomes here instead
  cd $PROJECTFOLDER/resources
  wget https://a3s.fi/Holoruminant-data/2024.09.18.reference.tar.gz
  tar -xvf 2024.09.18.reference.tar.gz

# For MAGScot are also dedicated files needed, which can be pulled in a similar way
  cd $PROJECTFOLDER/resources
  wget https://a3s.fi/Holoruminant-data/2024.09.18.MAGScot.tar.gz
  tar -xvf 2024.09.18.MAGScot.tar.gz

# Get the example read data
  cd $PROJECTFOLDER
  wget https://a3s.fi/Holoruminant-data/2024.09.18.reads.tar.gz
  tar -xvf 2024.09.18.reads.tar.gz
```

If you have downloaded the resources already into another project, you can share the resources also to a new project, e.g. by creating a symbolic link

```
cd $PROJECTFOLDER
ln -s /project/with/existing/resources resources

```

Now we copy the configuration files from the pipeline folder to the project folder, to adjust the configurations to the project specifics

```
cd $PIPELINEFOLDER
cp -r config $PROJECTFOLDER
cp -r run_Pipeline-Holoruminant-meta.sh $PROJECTFOLDER
```

# Setting up the pipeline

## run_Pipeline-Holoruminant-meta.sh
This is the pipeline starting wrapper script. It takes care of enabling Snakemake (e.g. in case you have it as a module on your server) and also wraps the Snakemake options nicely. Furthermore, it handles to setup the environment variables for tmp and cache folders of apptainer or singularity and also can be used to prepare the rulegraph.

Enter the required values and paths according to the comments in the file.

## config/config.yaml
Here are the paths to the different configuration files stored, which might not need any adjustments from the user (e.g. for Holoruminant users). 

In addition, the specs for the resource allocations are provided here. The defaults are currently not calibrated and need still some closer evaluation. Adjust the values to your needs and names from your hpc (like queue names)

## config/profiles/
Here are the HPC profiles stored. The current default configuration is adjusted to our system called Puhti and is located in the subfolder `Puhti/` in the file `config.yaml`. For your own system, create a new subfolder with the name of your system and copy the config file from `Puhti/` there to adjust. Please do not rename the yaml file, it needs to be `config.yaml`. 

Here, set your typical default resources and check what requriements your generic slurm (or whatever executor you use) command has. In essence, the `cluster-generic-submit-cmd` needs to match the requirements of your system, e.g. the `slurm_account` option might be very specific on our system Puhti and might not be accepted or required on your system.

## config/features.yaml
Here we can adjust the reference genomes and databases that should be used from the pipeline. The
current defaults are for the Holoruminant project and have as such a very specific set of reference genomes that are used for filtering and checking contaminations. Adjust yours in the `hosts:` section. Here, you can just use a own name, followed by collon and the path to it. Reference genomes are expected to be gzipped.

In the `databases:` section the paths to the corresponding databases are used. The default paths meet the folder structure you will obtain, when you download the databases from our server. The main adjustments to do are a) the kraken path and b) phylophlan. Here, sub-databases can be given and the tools run one after another the searches against these databases. In case of kraken, we provide a small standard database as well as a rather large one called `refseq500`. Just comment out the ones you do not want to use.

## config/params.yaml
This file contains the tuning parameters of the different tools. This file is far from being complete and is currently not calibrated. So, please check if the tools use the parameters you want them to use and add if needed the parameters to the rule and the params file.

## config/samples.tsv
This file contains the sample information and is required by the pipeline. 

The file has this structure

| sample_id  | library_id | forward_filename              | reverse_filename              | forward_adapter                           | reverse_adapter                           | assembly_ids |
|------------|------------|-------------------------------|-------------------------------|-------------------------------------------|-------------------------------------------|--------------|
| ERR2019410 | lib1       | reads/ERR2019410_1.fastq.gz   | reads/ERR2019410_2.fastq.gz   | AGATCGGAAGAGCACACGTCTGAACTCCAGTCA         | AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT         | ERR2019410   |
| ERR2019412 | lib1       | reads/ERR2019412_1.fastq.gz   | reads/ERR2019412_2.fastq.gz   | AGATCGGAAGAGCACACGTCTGAACTCCAGTCA         | AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT         | ERR2019412   |


where we can give the 

* `sample_id` , unique sample identified
* `library_id`, library identifiier for repeated libraries
* `forward_filename`, relative path to the forward read file of that sample
* `reverse_filename`, relative path to the reverse read file of that sample
* `forward_adapter`, adapter sequence from the forward read
* `reverse_adapter`, adapter sequence from the reverse read
* `assembly_ids`, group name for the co-assembly this sample should be used in

I am not sure, tbh, if currently several assembly_ids per sample are allowed, this would need to be tested.

There is a script in the script folder to create the sample sheet. For that, you can run

```
cd $PROJECTFOLDER
bash $PIPELINEFOLDER/workflow/scripts/createSampleSheet.sh
```

It should create the `samples.tsv` for the samples located in the `reads/` folder. You might need to adjust the script maybe accoring to the names of the reads or the adapter sequences you use.

In case you have several lanes for samples, you can concatenate them prior to creating the samples.tsv script with the script `concatenateFiles.sh`which is in the pipeline folder `workflow/scripts`. Currently, you would need to run the script inside the same folder where the fastq files are located.

# Usage
The pipeline can run the entire workflow at once. However, normally it is recommended to run different modules from the pipeline separated to get better control over the results and also to be able to react quicker to possible errors.

In the following it is assumed that the pipeline runs on a server that utilizes SLURM and Singularity. Other setups are also supported, but currently untested. In case you have a different setup and want to contribute a profile/configuration, please reach out.

For testing and developing, you can add to every command e.g. the option `-np` for a dry-run that prints the used commands.

The different module have also individual reports that can be generated by adding `report_` in front of the module name, when a module is called. However, the reports are currently under developments and do not produce any reasonable output and might crash even.

## 'reads-module'
Here some basic steps for the reads are performed.

Usage:

```
bash run_Pipeline-Holoruminant-meta.sh reads

# Generate the module report
bash run_Pipeline-Holoruminant-meta.sh report_reads
```

For testing, please check first the dry-run with commands printed by running the command like this

```
bash run_Pipeline-Holoruminant-meta.sh reads -np
```

For all other modules this works in a similar fashion, just add the `-np`-option for testing and `report_` to the module to generate a report for the module (not implemented yet for all modules).

## reference-module
The reference host genomes are recompressed in this module

Usage:

```
bash run_Pipeline-Holoruminant-meta.sh reference
```

## preprocess-module
This module covers the quality control and triming of the reads as well as the read-based
analysis and database searches.

Usage:

```
bash run_Pipeline-Holoruminant-meta.sh preprocess

# Generate the module report
bash run_Pipeline-Holoruminant-meta.sh report_preprocess
```

The command `preprocess` triggers the analysis of the entire preprocess submodule. However, individual tools can be called instead by running

```
bash run_Pipeline-Holoruminant-meta.sh preprocess__diamond
bash run_Pipeline-Holoruminant-meta.sh preprocess__krona
...
```

The underlying syntax is that an individual tool is always called by module name, followed by two underscores and the tool name. Please bear in mind that the pipeline automatically runs the required pre-steps. This syntax hold through-out the entire pipeline and can be used to run only subparts of the pipeline.


## assemble-module
This module runs all the assembly related tasks, like creating the metagenome and then the binning and combination of different binners.

Usage:

```
bash run_Pipeline-Holoruminant-meta.sh assemble
```

## quantify-module

After creating the assemblies, this module does the mapping of the reads and generates the quantification tables for the samples.

Usage:

```
bash run_Pipeline-Holoruminant-meta.sh quantify
```
## annotate-module

This module is the annotation workhorse, it aligns the mags against various databases and generates the annotation for them.

Usage:

```
bash run_Pipeline-Holoruminant-meta.sh annotate
```

## contig_annotate-module

This module is running the contig based analysis from beginning to end, meaning it first predicts genes on the metagenome assembly, annotates the predicted genes with eggnog, aligns the quality filtered reads to the assembly with bowtie and then quantifies it with coverm.

Usage:

```
bash run_Pipeline-Holoruminant-meta.sh contig_annotate
```

## report-module
Different modules create reports, here all reports at once can be generated. Otherwise, individual reports can also be reported after each module.

Usage:

```
bash run_Pipeline-Holoruminant-meta.sh report
```

# Module details

## 'reads' -module

The reads module can be started by

```
bash run_Pipeline-Holoruminant-meta.sh reads
```

and it triggers a set of rules:

### _reads__link
This rule makes a link to the original file, with a prettier name than default

It creates output files like this

```
forward_= results/reads/"{sample}.{library}_1.fq.gz",
reverse_= results/reads/"{sample}.{library}_2.fq.gz",
```

with `sample` and `library` being taken from the `sample.tsv` file.

### _helpers__fastqc
This rules creates fastqc reports for the input files and stores them also
in the results/reads folder as `*.zip` and `*.html` files.

## 'reference' - module
The reference module can be started by

```
bash run_Pipeline-Holoruminant-meta.sh reads
```

This module has essentially the only function to take the gezipped host genomes and
rezip them to bgzip

### _reference__hosts__recompress
This function loops through the provided host genomes and then recompresses them.

The results are then stored under `results/reference/hosts/`

## 'preprocess' - module
A larger subworkflow that consists of several steps. It can be run by

```
bash run_Pipeline-Holoruminant-meta.sh preprocess
```

### fastp subworflow
rule _preprocess__fastp__run:
    Run fastp on one library. THis is essentially the quality and adapter trimming. The output of quality trimming
    is then fed into the host decontamination.

### kraken2 subworkflow
TODO: CHECK; WHAT DATABASE TO USE?! https://github.com/R-Wright-1/kraken_metaphlan_comparison/wiki/Downloading-databases

wget https://genome-idx.s3.amazonaws.com/kraken/k2_nt_20231129.tar.gz

https://github.com/R-Wright-1/kraken_metaphlan_comparison/wiki/Downloading-databases

rule _preprocess__kraken2__assign:
    """
    Run kraken2 over all samples. Here, also the fastp reads are used, so we have not removed host contamination, when we assign kraken2 to the reads. Different databases used with kraken can be added into the configuration file `features.yaml`, intented after the kraken2 entry

### bowtie_2 subworkflow

This subworflow does a couple of things:

rule _preprocess__bowtie2__build:
    Build PRE_BOWTIE2 index for the host reference

rule _preprocess__bowtie2__map:
    Map one library to reference genome using bowtie2
    
rule _preprocess__bowtie2__extract_nonhost:
    Keep only pairs unmapped to the human reference genome, sort by name rather
    than by coordinate, and convert to FASTQ.

### Nonpareil subworkflow
rule _preprocess__nonpareil__run:
    Run nonpareil over one sample. For this step, the host decontaminated reads are used!
    Rodriguez-R LM & Konstantinidis KT (2014). Nonpareil: A redundancy-based approach to assess the level of coverage in metagenomic datasets. Bioinformatics 30 (5): 629-635.
    
### samtools subworflow
rule _preprocess__samtools__stats_cram:
    Compute the stats of a cram file using samtools stats. Here, we'll get mapping statistics for the host alignemtns.
    
### singlem subworkflow
TODO: CHECK WHAT DATABASE TO USE!!!

This part consists again out of several steps.. In details, these are

SingleM is a tool for profiling shotgun metagenomes. It has a particular strength in detecting microbial lineages which are not in reference databases. The method it uses also makes it suitable for some related tasks, such as assessing eukaryotic contamination, finding bias in genome recovery, computing ecological diversity metrics, and lineage-targeted MAG recovery.

Documentation can be found at https://wwood.github.io/singlem/

Citation

SingleM and Sandpiper: Robust microbial taxonomic profiles from metagenomic data. Ben J Woodcroft, Samuel T. N. Aroney, Rossen Zhao, Mitchell Cunningham, Joshua A. M. Mitchell, Linda Blackall, Gene W Tyson. bioRxiv 2024.01.30.578060; doi: https://doi.org/10.1101/2024.01.30.578060

rule _preprocess__singlem__pipe:
    Run singlem over one sample

rule _preprocess__singlem__condense:
    Aggregate all the singlem results into a single table
    
rule _preprocess__singlem__microbial_fraction:
    Run singlem microbial_fraction over one sample
    
rule _preprocess__singlem__aggregate_microbial_fraction:
    Aggregate all the microbial_fraction files into one tsv
    
### metaphlan4
For metaphlan 4 I downloaded the database like this

```
metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202403 --bowtie2db metaphlan4/
 metaphlan ERR2019411_1.fastq.gz,ERR2019411_2.fastq.gz --bowtie2out metagenomeERR.bowtie2.bz2 --nproc 20 --input_type fastq -o profiled_metagenome.txt --bowtie2db metaphlan4/
```
    
## 'assemble' -module

The assemble module runs first several assemblers and combines then the results. It can be initiated by running

```
bash run_Pipeline-Holoruminant-meta.sh assemble
```
    
It contains of several subworkflows, as

### megahit subworkflow
rule _assemble__megahit:
    Run megahit over samples associated to assembly, merging all libraries in the process. This creates the co-assemblies as defined in the samples.tsv file    

### bowtie2 - subworkflow

rule _assemble__bowtie2__build:
    Index a megahit assembly. Here, we prepare the megahit assembly for mapping
    
rule _assemble__bowtie2__map:
    Map one sample to one megahit assembly. Here, we map then the samples to the megahit assemblies
    
### concoct subworkflow

rule _assemble__concoct:
     This one takes all available assemblies and runs concot for them

### maxbin2 subworkflow
rule _assemble__maxbin2__run:
    """Run MaxBin2 over a single assembly, so for each assembly, we get the bins produced from maxbin2

### metabat2 subworkflow
rule _assemble__metabat2__run:
    Run metabat2 end-to-end on a single assembly, for all assemblies then


### drep subworflow
rule _assemble__drep__separate_bins:
    This one takes the magscot output and separates the bins
    
rule _assemble__drep__run:
    Dereplicate all the bins using dRep. This is the depreplication step for the bins

rule _assemble__drep__join_genomes:
    Join all the dereplicated genomes into a single file.

### magscot subworkflow
rule _assemble__magscot__prodigal:
    Run prodigal over a single assembly. So, we predict for each assembly the genes.    

rule _assemble__magscot__hmmsearch_pfam:
    Run hmmsearch over the predicted proteins of an assembly using Pfam as database  
  
rule _assemble__magscot__hmmsearch_tigr:
    Run hmmsearch over the predicted proteins of an assembly using TIGR as database  
CHECK HOW TO GET THE TIGR DATABASE!!!

rule _assemble__magscot__join_hmms:
    Join the results of hmmsearch over TIGR and Pfam
  
rule _assemble__magscot__merge_contig_to_bin:
    Merge the contig to bin files from CONCOCT, MaxBin2 and MetaBAT2  

rule _assemble__magscot__run:
    Run MAGSCOT over one assembly
  
rule _assemble__magscot__rename:
    Rename the contigs in the assembly to match the assembly and bin names  
  
### vamb subworkflow

## 'quantify' - module  
The quantify module quantifies the different assemblies and bins. You can use it like this

```
bash run_Pipeline-Holoruminant-meta.sh quantify
```

### bowtie2 subworkflow

rule _quantify__bowtie2__build:
    Index dereplicated genomes /mag
    
rule _quantify__bowtie2__map:
    Align one sample to the dereplicated genomes

### coverm subworkflow

rule _quantify__coverm__genome:
    Run CoverM genome for one library and one mag catalogue

rule _quantify__coverm__genome_aggregate:
    Run coverm genome and a single method

rule _quantify__coverm__contig:
    Run coverm contig for one library and one mag catalogue
    
rule _quantify__coverm__contig_aggregate:
    Run coverm contig and a single method
    
### samtools subworkflow
rule _quantify__samtools__stats_cram:
    Get stats from CRAM files using samtools stats.
    
## 'annotate' module    
This module takes care of the annotation steps

```
bash run_Pipeline-Holoruminant-meta.sh annotate
```
    
### quast subworkflow
rule annotate__quast:
    Run quast over one the dereplicated mags
    
### GTDBtk subworkflow
rule _annotate__gtdbtk__classify:
    Run GTDB-Tk over the dereplicated genomes

TODO: THE DATABASE NEEDS TO BE FETCHED AUTOMATICALLY, CURRENTLY THE USER NEEDS TO DO IT!    
TODO: UPDATE TO 2.4.0 (use release 220 instead)
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz
    
### DRAM subworkflow
rule _annotate__dram__annotate:
    Annotate dereplicate genomes with DRAM

rule _annotate__dram__distill:
    Distill DRAM annotations    

TODO: ALSO HERE THE DB NEEDS TO BE INSTALLED MANUALLY!
This needs quite much resources, so run it best in a sbatch script


#!/bin/bash
#SBATCH --job-name=create_dram_db
#SBATCH --account=project_2009831
#SBATCH --partition=hugemem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=800G
#SBATCH --time=24:00:00
#SBATCH --output=create_dram_db.out
#SBATCH --error=create_dram_db.err

singularity shell -B /scratch,/projappl,/users,/dev/shm:/tmp,/run:/run .snakemake/singularity/675d014754a2522c1e382e4b6f21b014.simg
DRAM-setup.py export_config > my_old_config.txt
export DRAM_CONFIG_LOCATION=my_old_config.txt
DRAM-setup.py prepare_databases --output_dir resources/databases/dram/20230811/

### CheckM2 subworkflow
rule _annotate__checkm2__predict:
    Run CheckM2 over the dereplicated mags
Install the database manually first

# DRAM database generation
Apptainer> DRAM-setup.py export_config > /scratch/project_2009831/Pipe_dev/my_dram_config.txt
Apptainer> export DRAM_CONFIG_LOCATION=/scratch/project_2009831/Pipe_dev/my_dram_config.txt
Apptainer> DRAM-setup.py prepare_databases --output_dir /scratch/project_2009831/Pipe_dev/resources/databases/dram/20230811/

https://github.com/WrightonLabCSU/DRAM/issues/26#issuecomment-685212290

# Bakta
The bakta database I downloaded is 5.1:
https://zenodo.org/records/10522951


## References

- [`fastp`](https://github.com/OpenGene/fastp)
- [`kraken2`](https://github.com/DerrickWood/kraken2)
- [`SingleM`](https://github.com/wwood/singlem)
- [`Nonpareil`](https://github.com/lmrodriguezr/nonpareil)
- [`bowtie2`](https://github.com/BenLangmead/bowtie2)
- [`samtools`](https://github.com/samtools/samtools)
- [`MEGAHIT`](https://github.com/voutcn/megahit)
- [`CONCOCT`](https://github.com/BinPro/CONCOCT)
- [`MaxBin2`](http://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html)
- [`MetaBat2`](https://bitbucket.org/berkeleylab/metabat)
- [`MAGScoT`](https://github.com/ikmb/MAGScoT)
- [`dRep`](https://github.com/MrOlm/drep)
- [`QUAST`](https://github.com/ablab/quast)
- [`GTDB-TK`](https://github.com/Ecogenomics/GTDBTk)
- [`DRAM`](https://github.com/WrightonLabCSU/DRAM)
- [`CoverM`](https://github.com/wwood/CoverM)
- [`FastQC`](https://github.com/s-andrews/FastQC)
- [`multiqc`](https://github.com/ewels/MultiQC)

# Acknowledgements
This pipeline is a fork from the Snakemake workflow

https://github.com/3d-omics/mg_assembly/

and tailored and extended to the needs of the Holoruminant project.
