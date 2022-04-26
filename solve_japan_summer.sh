#!/bin/bash

#SBATCH --job-name=Msummer
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpuonly
#SBATCH --mem=16G
#SBATCH --time=50000

echo "Starting job on node $SLURMD_NODENAME"
cd /users/mlevorato/doutorado/manutencao_pocos/src/julia
sh ./install_gurobi_license.sh $SLURMD_NODENAME

echo "CCP Julia Japan spring experiment script start."

cd /users/mlevorato/doutorado/robusto/RCCP
/users/mlevorato/julia-1.6.0/bin/julia solve_models_japan_summer.jl
