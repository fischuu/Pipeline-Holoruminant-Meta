# Set the project relevant paths
################################################################################
pipelineFolder="/users/fischerd/git/Pipeline-Holoruminant-Meta"
projectFolder="/scratch/project_2009831/Pipe_dev"

# Setup the server profile. You can copy the default profile and adjust it to your hpc.
# Check that the generic command meets the requirements from your (slurm) executor.
# It is recommended to copy the Profile folder to a new folder with the name of
# your cluster
################################################################################
Profile=$pipelineFolder/config/profiles/Puhti

# Make Snakemake available. Here is an example, how it would be loaded with the
# module system, but other systems might have it available per standard installation.
# In that case you can comment out the follwing line
################################################################################
module load snakemake/8.4.6

# For use with Apptainer, set these variables
################################################################################
export APPTAINER_TMPDIR="/scratch/project_2009831/tmp"
export APPTAINER_CACHEDIR="/scratch/project_2009831/tmp"
mkdir -p $APPTAINER_TMPDIR
mkdir -p $APPTAINER_CACHEDIR

# For use with Singularity, set these variables
################################################################################
# export SINGULARITY_TMPDIR="/scratch/project_2009831/tmp"
# export SINGULARITY_CACHEDIR="/scratch/project_2009831/tmp"
# mkdir -p $SINGULARITY_TMPDIR
# mkdir -p $SINGULARITY_CACHEDIR

# Create the rulegraph
################################################################################
#snakemake -s $pipelineFolder/workflow/Snakefile \
#          --configfile $projectFolder/config/config.yaml \
#          --rulegraph | dot -T png > $projectFolder/workflow.png


# Run the pipeline (singularity/apptainer)
# Important adjustments: 
# --jobs, how many parallel jobs do you allow on your system for this project
# --singularity-args, are the important folders and paths bound to the containers? If you are unsure, go with the defaults and check for errors
################################################################################
snakemake -s $pipelineFolder/workflow/Snakefile \
          --jobs 150 \
          --use-singularity \
          --configfile $projectFolder/config/config.yaml \
          --profile $Profile \
          --singularity-args "-B /scratch,/projappl,/users,/dev/shm:/tmp,/run" \
          --latency-wait 60 \
          --scheduler greedy \
          $@
