# Generate more configuration keys
# Define the script_folder dynamically based on the pipeline_folder
SCRIPT_FOLDER = os.path.join(config["pipeline_folder"], "workflow", "scripts")


READS = Path("results/reads/")
WD = os.getcwd()

# reference
REFERENCE = Path("results/reference/")
HOSTS = REFERENCE / "hosts"

# preprocess
PRE = Path("results/preprocess/")
FASTP = PRE / "fastp/"
PRE_INDEX = PRE / "index/"
PRE_BOWTIE2 = PRE / "bowtie2"
NONHOST = PRE / "nonhost/"
PRE_QUANT = PRE / "quantification/"

# preprocess evaluation
NONPAREIL = PRE / "nonpareil/"
METAPHLAN = PRE / "metaphlan/"
HUMANN = PRE / "humann/"
PHYLOFLASH = PRE / "phyloflash/"
SINGLEM = PRE / "singlem/"
PRE_COVERM = PRE / "coverm/"
KRAKEN2 = PRE / "kraken2/"
KRONA = PRE / "krona/"
DIAMOND = PRE / "diamond/"

# assemble
ASSEMBLE = Path("results/assemble/")
MEGAHIT = ASSEMBLE / "megahit/"
METASPADES = ASSEMBLE / "metaspades/"
ASSEMBLE_RENAME = ASSEMBLE / "renaming/"
ASSEMBLE_INDEX = ASSEMBLE / "index/"
ASSEMBLE_BOWTIE2 = ASSEMBLE / "bowtie2/"
CONCOCT = ASSEMBLE / "concoct/"
METABAT2 = ASSEMBLE / "metabat2/"
MAXBIN2 = ASSEMBLE / "maxbin2/"
# VAMB = ASSEMBLE / "vamb/"  # This could be a metabinner
MAGSCOT = ASSEMBLE / "magscot/"
PRODIGAL = MAGSCOT / "prodigal/"
DREP = ASSEMBLE / "drep/"

# quantify
QUANT = Path("results/quantify/")
QUANT_INDEX = QUANT / "index/"
QUANT_BOWTIE2 = QUANT / "bowtie2/"
COVERM = QUANT / "coverm/"

# dereplicate evaluation
ANN = Path("results/annotate/")
GTDBTK = ANN / "gtdbtk/"
QUAST = ANN / "quast/"
CAMPER = ANN / "camper/"
DRAM = ANN / "dram/"
CHECKM = ANN / "checkm2"
BAKTA = ANN / "bakta"
EGGNOG = ANN / "eggnog"
PHYLOPHLAN = ANN / "phylophlan"
PROTEINORTHO = ANN / "proteinortho"
SYLPH = ANN / "sylph"

# contig annotation
CONTIG = Path("results/contig_annotate/")
CONTIG_PRODIGAL = CONTIG / "prodigal"
CONTIG_EGGNOG = CONTIG / "eggnog"
CONTIG_FEATURECOUNTS = CONTIG / "featurecounts"

# reports
REPORT = Path("reports/")
PIPELINE_REPORT = Path("Rreports/")
REPORT_STEP = REPORT / "step/"
REPORT_SAMPLE = REPORT / "library/"
