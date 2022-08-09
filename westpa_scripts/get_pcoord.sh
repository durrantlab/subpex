#!/bin/bash

# This script is run when calculating initial progress coordinates for new
# initial states (istates).  This script is NOT run for calculating the progress
# coordinates of most trajectory segments; that is instead the job of runseg.sh.

export ENGINE="NAMD"

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

# Link simulation files
[[ ! -e mol.pdb ]] && ln -s $WEST_SIM_ROOT/reference/mol.pdb mol.pdb
[[ ! -e mol.prmtop ]] && ln -s $WEST_SIM_ROOT/reference/mol.prmtop mol.prmtop
if [ "$ENGINE" == "NAMD" ] ; then
  [[ ! -e seg.dcd ]] && ln -s $WEST_SIM_ROOT/reference/mol.dcd seg.dcd
elif [ "$ENGINE" == "AMBER" ] ; then
  [[ ! -e seg.rst ]] && ln -s $WEST_SIM_ROOT/reference/mol.nc seg.rst
fi

# Use a custom script to calculate the jaccard distance between the starting
# structure and the initial state (should be 0 since we are copying the files).

python3 $WEST_SIM_ROOT/westpa_scripts/pcoord_istate.py mol.pdb $WEST_SIM_ROOT/west.cfg

cp pcoord.txt $WEST_PCOORD_RETURN  # Always exists
[[ -e pvol.txt ]] && cp pvol.txt $WEST_PVOL_RETURN
[[ -e rog.txt ]] && cp rog.txt $WEST_ROG_RETURN
[[ -e bb.txt ]] && cp bb.txt $WEST_BB_RETURN
[[ -e composite.txt ]] && cp composite.txt $WEST_COMPOSITE_RETURN
[[ -e prmsd.txt ]] && cp prmsd.txt $WEST_PRMSD_RETURN
[[ -e jd.txt ]] && cp jd.txt $WEST_JD_RETURN

# If we are running in debug mode, then output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
