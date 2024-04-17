#!/bin/bash
############################## Submit Job in Julia ######################################
#SBATCH --time=48:00:00
#SBATCH --job-name="bcd"
#SBATCH --mail-user=yalan@stanford.edu
#SBATCH --mail-type=END
#SBATCH --output=bcd_e%j.txt
#SBATCH --error=FAILURE_bcd_e%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --partition=normal
#####################################

# Load module for Gurobi and Julia (should be most up-to-date version, i.e. 1.7.2)
module load julia
module load gurobi

# Change to the directory of script
# export SLURM_SUBMIT_DIR=/home/users/yalan/binary_seq_opt/sherlock_demo
export SLURM_SUBMIT_DIR=/home/groups/gracegao/prn_codes/binary_seq_opt/sherlock_demo

# Change to the job directory
cd $SLURM_SUBMIT_DIR

export GUROBI_HOME="/share/software/user/restricted/gurobi/10.0.1_py39"
# export MOSEKBINDIR="/home/groups/gracegao/mosek/mosek/9.3/tools/platform/linux64x86/bin"

lscpu

# Run script
# julia --heap-size-hint=4G sherlock_demo.jl 0 127 66 20 SOS 10 100000 100000 10 8 false
julia --heap-size-hint=8G sherlock_demo.jl 0 1023 31 15 MSOS false 5 100000 100000 10 2 false 3 true
julia --heap-size-hint=8G sherlock_demo.jl 0 1023 31 15 SOS false 5 100000 100000 10 2 false 3 true
