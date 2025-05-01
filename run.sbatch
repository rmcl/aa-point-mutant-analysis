#!/bin/bash
#SBATCH --job-name=mutation_analysis
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err
#SBATCH --array=0-40
#SBATCH --partition=memory
#SBATCH --ntasks=1                 # Number of tasks per node
#SBATCH --cpus-per-task=1          # Number of CPU cores per task
#SBATCH --time=12:00:00            # Time limit
#SBATCH --mem=8G                   # Memory per node

# Load necessary modules (if applicable)
# module load rosetta

OUTPUT_DIR=outputs/

source /packages/miniconda3/20240410/etc/profile.d/conda.sh
conda activate rosetta_env


# Set up the environment
export LD_LIBRARY_PATH=./:$LD_LIBRARY_PATH

# Calculate the partition
pdb_list="inputs.txt"
total_files=$(wc -l < $pdb_list)
files_per_task=$(( (total_files + SLURM_ARRAY_TASK_COUNT - 1) / SLURM_ARRAY_TASK_COUNT ))  # Round up
start_index=$(( SLURM_ARRAY_TASK_ID * files_per_task ))
end_index=$(( start_index + files_per_task - 1 ))

# Process the assigned files
echo "Task ID: $SLURM_ARRAY_TASK_ID, Processing files $start_index to $end_index"
sed -n "$((start_index + 1)),$((end_index + 1))p" $pdb_list | while read mutation_detail; do
    echo "MODIFY MUTATION: ${mutation_detail}"

    python ./single_aa_score_batch.py ${mutation_detail}

done
