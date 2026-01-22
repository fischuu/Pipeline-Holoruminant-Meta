# Troubleshooting

## Starting the pipeline

_The reads-module does not even finish correctly_

1. Check the slurm output in .snakemake/slurm_logs, can you find something there?
2. Check the binding folders (-B) in run_Pipeline-Holoruminant-meta.sh, do they all exist?