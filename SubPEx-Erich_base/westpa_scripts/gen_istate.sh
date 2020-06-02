#!/bin/bash


# For debugging add to the script --debug
if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

# Custom script to calculate the reference
python3 $WEST_SIM_ROOT/westpa_scripts/get_reference_fop.py settings.json

#Make the directories for each of the initial states and symlink the necessary files.
cd $WEST_SIM_ROOT

mkdir -p $(dirname $WEST_ISTATE_DATA_REF)
ln -s $WEST_BSTATE_DATA_REF $WEST_ISTATE_DATA_REF

