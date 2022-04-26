#!/bin/bash

#SBATCH --job-name=CCPJwinter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpuonly
#SBATCH --mem=32G
#SBATCH --time=50000

echo "Starting job on node $SLURMD_NODENAME"
cd /users/mlevorato/doutorado/manutencao_pocos/src/julia
sh ./install_gurobi_license.sh $SLURMD_NODENAME

echo "CCP Julia Japan winter experiment script start."

cd /users/mlevorato/doutorado/robusto/RCCP
/users/mlevorato/julia-1.6.0/bin/julia --threads=8 run_sim_instance_deltamin10_winter.txt.jl
