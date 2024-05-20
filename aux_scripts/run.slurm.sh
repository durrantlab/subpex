#!/bin/bash
#SBATCH --job-name=subpex_1
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --output=job_logs/slurm.out
#SBATCH --error=job_logs/slurm.err
#SBATCH --time=24:00:00
#SBATCH --cluster=mpi
#SBATCH --partition=opa-high-mem
#SBATCH --mail-user=user@email.domain
#SBATCH --mail-type=END,FAIL  #BEGIN

# TODO: Erich, can we get example slurm scripts for AMBER on GPU too?
# TODO: Look to GPU scripts in env.sh (./env-crc/) for inspiration. Not AMBER/NAMD specific.

#
# Example script for submitting a weighted ensemble simulation to the slurm job
# scheduler. Be sure to modify the directives above and to run init.sh first! If
# you want to run SubPEx directly from the command line, use run.sh.
#

source env.sh

SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info-$SLURM_JOBID.json
echo $WEST_PYTHON

w_run --work-manager=zmq --n-workers=0 --zmq-mode=master --zmq-write-host-info=$SERVER_INFO --zmq-comm-mode=tcp &> ./job_logs/west-$SLURM_JOBID.log &

# wait on host info file up to one minute
for ((n=0; n<60; n++)); do
    if [ -e $SERVER_INFO ] ; then
        echo "== server info file $SERVER_INFO =="
        cat $SERVER_INFO
        break
    fi
    sleep 1
done

# exit if host info file doesn't appear in one minute
if ! [ -e $SERVER_INFO ] ; then
    echo 'server failed to start'
    exit 1
fi

# start clients, with the proper number of cores on each

scontrol show hostname $SLURM_NODELIST > slurm_nodelist.txt

for node in $(scontrol show hostname $SLURM_NODELIST); do
    ssh -o StrictHostKeyChecking=no $node $PWD/node.sh $SLURM_SUBMIT_DIR $SLURM_JOBID $node $CUDA_VISIBLE_DEVICES --work-manager=zmq --n-workers=1 --zmq-mode=client --zmq-read-host-info=$SERVER_INFO --zmq-comm-mode=tcp &
done

wait
