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

# If using NAMD...
ln -s $WEST_SIM_ROOT/reference/{REFERENCE}.pdb mol.pdb
ln -s $WEST_SIM_ROOT/reference/{RESTART_FILE} seg.dcd  # TODO: Should this be an rst file?
ln -s $WEST_SIM_ROOT/reference/{TOPOLOGY_FILE} mol.prmtop

# If using AMBER... TODO: Can Erich fill in below?
# ln -s $WEST_SIM_ROOT/reference/{REFERENCE}.pdb mol.pdb
# ln -s $WEST_SIM_ROOT/reference/{RESTART_FILE} seg.dcd
# ln -s $WEST_SIM_ROOT/reference/{TOPOLOGY_FILE} mol.prmtop

# Use a custom script to calculate the jaccard distance between the starting
# structure and the initial state (should be 0 since we are copying the files).

python3 $WEST_SIM_ROOT/westpa_scripts/pcoord_istate.py mol.pdb $WEST_SIM_ROOT/west.cfg

# TODO: Erich, can we just always copy over all possible aux data, regardless of
# what specified? That way we don't have to ask users to edit this portion. Let
# me know if that seems reasonable.
[[ -e pcoord.txt ]] && cp pcoord.txt $WEST_PCOORD_RETURN
[[ -e pvol.txt ]] && cp pvol.txt $WEST_PVOL_RETURN
[[ -e rog.txt ]] && cp rog.txt $WEST_ROG_RETURN
[[ -e bb.txt ]] && cp bb.txt $WEST_BB_RETURN
[[ -e fop.txt ]] && cp fop.txt $WEST_FOP_RETURN
# TODO: prmsd? rog_pocket instead of rog? jd? Any other possible auxdata
# possible here?

# If we are running in debug mode, then output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
