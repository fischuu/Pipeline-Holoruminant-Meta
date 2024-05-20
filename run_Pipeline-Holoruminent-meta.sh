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

# For use of Apptainer, export these variables
################################################################################
#export SINGULARITY_TMPDIR="/scratch/project_2001746/tmp"
#export SINGULARITY_CACHEDIR="/scratch/project_2001746/tmp"
#mkdir -p $SINGULARITY_TMPDIR

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

#          --cluster-config $projectFolder/Snakebite-GBS_server-config.yaml \
#          --cluster "sbatch -t {cluster.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type} -p {cluster.partition}" \
