#!/bin/bash

# This is a script to prepare the restart of a SubPEx simulation (any WESTPA sim should also work) 
# without loosing any information.

while getopts ":n:e:h" flag
do
    case "${flag}" in
        n)
            GEN=$OPTARG
            ;;
        e)
            conda activate $OPTARG
            ;;
        h)
            echo "Usage: restart_subpex.sh [options] -s [GENERATION]"
            echo " "
            echo "This is a scrpit to automate the truncation in the west.h5 it has to be run from WEST_HOME"
            echo "It does the following:"
            echo " - saving a copy of the original west.h5 for save keeping"
            echo " - moving the iterations starting from the trucated one to reference directory"
            echo " - moving the binbounds.txt file to a new one"
            exit 1
            ;;
    esac
done

cp west.h5 west_old.h5
w_truncate -n $GEN
mv binbounds.txt binbounds_old.txt

# Getting the las iteration that was run
LASTITER="$(ls -Art traj_segs | tail -n 1)"
LASTITERNUM=$((10#$LASTITER))

# Move all the gnerations that where trucated
for i in $(seq -f "%06g" $GEN $LASTITERNUM)
do 
    mv traj_segs/$i reference/
done
