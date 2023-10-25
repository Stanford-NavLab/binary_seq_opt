#!/bin/bash
############################## Submit Job in Julia ######################################
#SBATCH --time=48:00:00
#SBATCH --job-name="257timing_eurasip"
#SBATCH --mail-user=yalan@stanford.edu
#SBATCH --mail-type=END
#SBATCH --output=257timing_eurasip_e%j.txt
#SBATCH --error=FAILURE_257timing_eurasip_e%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH --partition=normal
#####################################

# Load module for Gurobi and Julia (should be most up-to-date version, i.e. 1.7.2)
module load julia
module load gurobi

# Change to the directory of script
export SLURM_SUBMIT_DIR=/home/users/yalan/binary_seq_opt/eurasip

# Change to the job directory
cd $SLURM_SUBMIT_DIR

# Run script
julia --heap-size-hint=8G timing.jl 0 257 130 24 SOS 10 1 1 50
julia --heap-size-hint=8G timing.jl 0 127 66 24 SOS 10 1 1 50
