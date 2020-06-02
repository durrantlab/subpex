#!/bin/bash
#SBATCH --job-name=subpex_1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output=job_logs/slurm.out
#SBATCH --error=job_logs/slurm.err
#SBATCH --time=128:00:00
#SBATCH --cluster=smp
#
# run.sh
#
# Run the weighted ensemble simulation. Make sure you ran init.sh first!
#

module purge
module load namd/2.13-multicore
module unload python
source env.sh
export WEST_PYTHON="/ihome/jdurrant/erh91/miniconda3/envs/py2_md/bin/python2.7"
echo $WEST_PYTHON 
#source env.sh

rm -f west.log
$WEST_ROOT/bin/w_run --work-manager processes "$@" &> west.log

#echo $WEST_PYTHON
