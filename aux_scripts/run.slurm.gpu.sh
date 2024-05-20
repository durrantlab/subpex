#!/bin/bash
#SBATCH --job-name=subpex_1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --output=job_logs/slurm.out
#SBATCH --error=job_logs/slurm.err
#SBATCH --time=24:00:00
#SBATCH --cluster=gpu
#SBATCH --partition=titanx
#SBATCH --mail-user=user@email.domain
#SBATCH --mail-type=END,FAIL  #BEGIN
#SBATCH --gres=gpu:1

# TODO: Erich, can we get example slurm scripts for AMBER on GPU too?
# TODO: Look to GPU scripts in env.sh for inspiration. Not AMBER/NAMD specific.

#
# Example script for submitting a weighted ensemble simulation to the slurm job
# scheduler. Be sure to modify the directives above and to run init.sh first! If
# you want to run SubPEx directly from the command line, use run.sh.
#

source env.sh

SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info-$SLURM_JOBID.json
echo $WEST_PYTHON

#w_run --work-manager=serial --n-workers=0 &> ./job_logs/west-$SLURM_JOBID.log
w_run --work-manager=serial &> ./job_logs/west-$SLURM_JOBID.log


# start clients, with the proper number of cores on each

scontrol show hostname $SLURM_NODELIST > slurm_nodelist.txt
