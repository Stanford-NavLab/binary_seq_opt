#!/bin/bash
############################## Submit Job in Julia ######################################
#SBATCH --time=24:00:00
#SBATCH --job-name="bcd"
#SBATCH --mail-user=yalan@stanford.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=bcd_e%j.txt
#SBATCH --error=FAILURE_bcd_e%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH --partition=normal
#####################################

# Load module for Gurobi and Julia (should be most up-to-date version, i.e. 1.7.2)
module load julia
module load gurobi

# Change to the directory of script
export SLURM_SUBMIT_DIR=/home/users/yalan/binary_seq_opt/sherlock_demo

# Change to the job directory
cd $SLURM_SUBMIT_DIR

export GUROBI_HOME="/share/software/user/restricted/gurobi/9.0.3_py36"
# export MOSEKBINDIR="/home/groups/gracegao/mosek/mosek/9.3/tools/platform/linux64x86/bin"

lscpu

# Run script
julia sherlock_demo.jl 0 31 7 15 ISL
