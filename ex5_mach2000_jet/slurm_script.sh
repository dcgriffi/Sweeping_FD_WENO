#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 30:00:00
#SBATCH -o output/log_%j.out
#SBATCH -e output/log_%j.err

# Load modules (module names are site-specific; adjust for your cluster)
module purge
module load gcc/13.1.0

#Choose mode: testing or final
run_mode="testing"

# Set base directory names
output_root="${OUTPUT_ROOT:-./runs}"
timestamp=$(date +%Y%m%d_%H%M%S)
run_id="run_$timestamp"


if [ "$run_mode" == "testing" ]; then
    base_storage_dir="$output_root/testing_output/"
    else
    base_storage_dir="$output_root/final_output/"
    fi

# Ensure actual storage location exists
mkdir -p "$base_storage_dir"

# Make run-specific subfolder
run_output_dir="$base_storage_dir/$run_id"
mkdir -p "$run_output_dir"

# Run the executable, passing the output directory path
# (this serial code does not read the argument; it writes its .dat files to the working directory)
srun ./a "$run_output_dir"
