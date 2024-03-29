#!/bin/bash

# This script prepares the restart of a SubPEx simulation (or any WESTPA
# simulation) without loosing any information.

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
            echo "Usage: restart_subpex.sh [options] -n [GENERATION]"
            echo " "
            echo "This script automates the truncation of the west.h5 file. It must be run from WEST_HOME"
            echo "It does the following:"
            echo " - saves a copy of the original west.h5 for safe keeping"
            echo " - moves the iterations starting from the trucated one to the reference directory"
            echo " - moves the binbounds.txt file to a new file"
            exit 1
            ;;
    esac
done

# If west.h5 exists
if [ -f west.h5 ]; then
    cp west.h5 west_old.h5
fi
w_truncate -n $GEN

if [ -f binbounds.txt ]; then
    mv binbounds.txt binbounds_old.txt
fi

# Getting the last iteration that finished
LASTITER="$(ls -Art traj_segs | tail -n 1)"
LASTITERNUM=$((10#$LASTITER))

# Move all the generations that were trucated
for i in $(seq -f "%06g" $GEN $LASTITERNUM)
do
    TRAJSEGDEST="reference/$i"
    # Does the directory exist?
    if [ -d $TRAJSEGDEST ]; then
        # Move the file, .BAK
        mv "$TRAJSEGDEST" "$TRAJSEGDEST".BAK.$(date +%s)
    fi
    mv traj_segs/$i $TRAJSEGDEST
done
