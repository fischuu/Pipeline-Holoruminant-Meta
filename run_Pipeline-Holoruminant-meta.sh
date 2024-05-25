# Set the project relevant paths
################################################################################
pipelineFolder="/users/fischerd/git/Pipeline-Holoruminant-Meta"
projectFolder="/scratch/project_2009831/Pipe_dev"
Profile=$pipelineFolder/config/profiles/Puhti

# Make Snakemake available 
################################################################################
module load snakemake/8.4.6

# For use with Apptainer, export these variables
################################################################################
export APPTAINER_TMPDIR="/scratch/project_2009831/tmp"
export APPTAINER_CACHEDIR="/scratch/project_2009831/tmp"
mkdir -p $APPTAINER_TMPDIR

export SINGULARITY_TMPDIR="/scratch/project_2009831/tmp"
export SINGULARITY_CACHEDIR="/scratch/project_2009831/tmp"

# Create the rulegraph
################################################################################
#snakemake -s $pipelineFolder/workflow/Snakefile \
#          --configfile $projectFolder/Snakebite-GBS_config.yaml \
#          --rulegraph | dot -T png > $projectFolder/workflow.png


# Run the pipeline (singularity)
################################################################################
snakemake -s $pipelineFolder/workflow/Snakefile \
          --jobs 150 \
          --use-singularity \
          --configfile $projectFolder/config/config.yaml \
          --profile $Profile \
          --singularity-args "-B /scratch,/projappl,/users,/dev/shm:/tmp,/run:/run" \
          --latency-wait 60 \
          --scheduler greedy \
          $@
