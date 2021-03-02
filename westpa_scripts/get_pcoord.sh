#!/bin/bash

# This script is run when calculating initial progress coordinates for new
# initial states (istates).  This script is NOT run for calculating the progress
# coordinates of most trajectory segments; that is instead the job of runseg.sh.

# If we are debugging, output a lot of extra information. Option --debug
if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

# Run the script to obtain hte refernece fop and selection string.
#cd $WEST_SIM_ROOT
#python $WEST_SIM_ROOT/westpa_scripts/get_reference_fop.py $WEST_SIM_ROOT/reference/settings.yaml
#cd $WEST_SIM_ROOT

# Make sure we are in the correct directory
cd $WEST_SIM_ROOT
source env.sh
cd $WEST_STRUCT_DATA_REF

ln -s $WEST_SIM_ROOT/reference/mol.pdb .
ln -s $WEST_SIM_ROOT/reference/mol.prmtop . 
ln -s $WEST_SIM_ROOT/reference/mol.inpcrd .

ln -s $WEST_SIM_ROOT/reference/seg.coor .
ln -s $WEST_SIM_ROOT/reference/seg.dcd .
ln -s $WEST_SIM_ROOT/reference/seg.vel .
ln -s $WEST_SIM_ROOT/reference/seg.xsc .

#ln -s $WEST_SIM_ROOT/reference/ref.pdb mol.pdb
#ln -s $WEST_SIM_ROOT/reference/equil_npt.rst seg.rst
#ln -s $WEST_SIM_ROOT/reference/mol.prmtop . 


# Use a custom script to calculate the jaccard distance between the starting
# structure and the initial state (should be 0 since we are copying the files).

python3 $WEST_SIM_ROOT/westpa_scripts/pcoord_istate.py mol.pdb $WEST_SIM_ROOT/reference/settings.cfg

cp pcoord.txt $WEST_PCOORD_RETURN
cp pvol.txt $WEST_PVOL_RETURN
cp rog.txt $WEST_ROG_RETURN
cp bb_rmsd.txt $WEST_BB_RETURN
cp fop.txt $WEST_FOP_RETURN

# If we are running in debug mode, then output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
