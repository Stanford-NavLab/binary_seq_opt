#!/bin/bash
############################## Submit Job in Julia ######################################
#SBATCH --time=48:00:00
#SBATCH --job-name="bcd_icassp23"
#SBATCH --mail-user=yalan@stanford.edu
#SBATCH --mail-type=END
#SBATCH --output=bcd_icassp23_e%j.txt
#SBATCH --error=FAILURE_bcd_icassp23_e%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=2G
#SBATCH --partition=normal
#####################################

# Load module for Gurobi and Julia (should be most up-to-date version, i.e. 1.7.2)
module load julia
module load gurobi

# Change to the directory of script
export SLURM_SUBMIT_DIR=/home/users/yalan/binary_seq_opt/icassp23

# Change to the job directory
cd $SLURM_SUBMIT_DIR

export GUROBI_HOME="/share/software/user/restricted/gurobi/9.0.3_py36"
# export MOSEKBINDIR="/home/groups/gracegao/mosek/mosek/9.3/tools/platform/linux64x86/bin"

lscpu

# Run script
julia icassp23.jl 0 ../results/BCD-SOS-BiST_64_4-14_34_57_3457.jls 64 4 20 SOS false 256 256 false
# julia icassp23.jl 0 ../results/BCD-SOS-BiST_64_4-14_35_14_3514.jls 64 4 20 SOS false 256 256 false
# julia icassp23.jl 0 ../results/BCD-SOS-BiST_64_4-14_35_33_3533.jls 64 4 20 SOS false 256 256 false
# julia icassp23.jl 0 ../results/BCD-SOS-BiST_64_4-14_35_50_3550.jls 64 4 20 SOS false 256 256 false
# julia icassp23.jl 0 ../results/BCD-SOS-BiST_64_4-14_36_08_368.jls 64 4 20 SOS false 256 256 false
# julia icassp23.jl 0 ../results/BCD-SOS-BiST_64_4-14_36_26_3626.jls 64 4 20 SOS false 256 256 false
# julia icassp23.jl 0 ../results/BCD-SOS-BiST_64_4-14_36_43_3643.jls 64 4 20 SOS false 256 256 false
# julia icassp23.jl 0 ../results/BCD-SOS-BiST_64_4-14_37_01_371.jls 64 4 20 SOS false 256 256 false
# julia icassp23.jl 0 ../results/BCD-SOS-BiST_64_4-14_37_18_3718.jls 64 4 20 SOS false 256 256 false
# julia icassp23.jl 0 ../results/BCD-SOS-BiST_64_4-14_37_36_3736.jls 64 4 20 SOS false 256 256 false
