#!/bin/bash
############################## Submit Job in Julia ######################################
#SBATCH --time=48:00:00
#SBATCH --job-name="127fix_eurasip"
#SBATCH --mail-user=yalan@stanford.edu
#SBATCH --mail-type=END
#SBATCH --output=127fix_eurasip_e%j.txt
#SBATCH --error=FAILURE_127fix_eurasip_e%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH --partition=normal
#####################################

# Load module for Gurobi and Julia (should be most up-to-date version, i.e. 1.7.2)
module load julia
module load gurobi

# Change to the directory of script
export SLURM_SUBMIT_DIR=/home/groups/gracegao/prn_codes/binary_seq_opt/eurasip

# Change to the job directory
cd $SLURM_SUBMIT_DIR

export GUROBI_HOME="/share/software/user/restricted/gurobi/9.0.3_py36"
# export MOSEKBINDIR="/home/groups/gracegao/mosek/mosek/9.3/tools/platform/linux64x86/bin"

lscpu

# Run script
julia --heap-size-hint=8G eurasip.jl 0 "" 127 66 25 ACZSOS false 1000000 1000000 false 66 1 false 100
