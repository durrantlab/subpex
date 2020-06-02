#!/bin/bash

# This script is run when calculating initial progress coordinates for new
# initial states (istates).  This script is NOT run for calculating the progress
# coordinates of most trajectory segments; that is instead the job of runseg.sh.

# If we are debugging, output a lot of extra information. Option --debug
if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

# Make sure we are in the correct directory
cd $WEST_SIM_ROOT
source env.sh
cd $WEST_STRUCT_DATA_REF

rm pcoord.txt pvol.txtbb_rmsd.txt rog.txt

ln $WEST_SIM_ROOT/reference/mol.pdb ref.pdb
ln $WEST_SIM_ROOT/reference/mol.pdb mol.pdb
ln $WEST_SIM_ROOT/westpa_scripts/settings.json .

ln -s $WEST_SIM_ROOT/reference/ref.coor .
ln -s $WEST_SIM_ROOT/reference/ref.dcd .
ln -s $WEST_SIM_ROOT/reference/ref.vel .
ln -s $WEST_SIM_ROOT/reference/ref.xsc .

# Use a custom script to calculate the jaccard distance between the starting
# structure and the initial state (should be 0 since we are copying the files).
python3 $WEST_SIM_ROOT/westpa_scripts/pcoord_istate.py ref.pdb mol.pdb settings.json --pvol --rog --bb_rmsd

cp pcoord.txt $WEST_PCOORD_RETURN
cp pvol.txt $WEST_PVOL_RETURN
cp rog.txt $WEST_ROG_RETURN
cp bb_rmsd.txt $WEST_BB_RETURN
cp fop.txt $WEST_FOP_RETURN

rm pcoord.txt pvol.txt rog.txt bb_rmsd.txt

# If we are running in debug mode, then output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
