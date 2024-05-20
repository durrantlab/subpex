#!/usr/bin/env bash

# This script prepares the restart of a SubPEx simulation (or any WESTPA
# simulation) without loosing any information.

# Assumes -n is always passed if script is called.

GEN=$1  # Get the generation number directly from the script argument

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
