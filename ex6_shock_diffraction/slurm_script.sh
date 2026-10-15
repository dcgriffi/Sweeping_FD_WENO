#!/bin/bash
#SBATCH -J MPI-bstep        # Job name
#SBATCH -n 32                       # Total number of MPI tasks (cut_x * cut_y = 32)
#SBATCH -N 1                        # Number of nodes (adjust if needed)
#SBATCH -t 17:00:00                 # Walltime (hh:mm:ss)
#SBATCH -o output/log_%j.out        # Standard output (create ./output before submitting)
#SBATCH -e output/log_%j.err        # Standard error

# -------------------------------
# Load environment
# (module names are site-specific; adjust for your cluster)
# -------------------------------

module purge
module load gcc/13.1.0-nvrtbp3     
module load hpcx-mpi/4.1.5rc2s-yflad4v

# -------------------------------
# User settings
# -------------------------------
run_mode="testing"   # Choose: "testing" or "final"

output_root="${OUTPUT_ROOT:-./runs}"

timestamp=$(date +%Y%m%d_%H%M%S)
run_id="run_$timestamp"

# Storage locations depend on mode
if [ "$run_mode" == "testing" ]; then
    base_storage_dir="$output_root/testing_output"
else
    base_storage_dir="$output_root/final_output"
fi

# Ensure output directories exist
mkdir -p "$base_storage_dir"
run_output_dir="$base_storage_dir/$run_id"
mkdir -p "$run_output_dir"

# -------------------------------
# Run MPI program
# -------------------------------
echo "Running job with $SLURM_NTASKS MPI ranks"
echo "Output directory: $run_output_dir"

srun ./a "$run_output_dir"
