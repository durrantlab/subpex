#!/bin/bash
#SBATCH --job-name=westpa-init
#SBATCH --output=westpa_init.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=12:00:00
#SBATCH --cluster=smp
#SBATCH --partition=high-mem

module load intel/2017.1.132
module load NAMD_MODULE
module load westpa/2017.10 

./init.sh
