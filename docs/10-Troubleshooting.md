# Troubleshooting

## Starting the pipeline

*The reads-module does not even finish correctly*

1. Check the slurm output in .snakemake/slurm_logs, can you find something there?
2. Check the binding folders (-B) in run_Snakebite-Holoruminant-MetaG.sh, do they all exist?