#!/bin/bash
#
# get_pcoord.sh
#
# This script is run when calculating initial progress coordinates for new 
# initial states (istates).  This script is NOT run for calculating the progress
# coordinates of most trajectory segments; that is instead the job of runseg.sh.

# If we are debugging, output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

# Make sure we are in the correct directory
cd $WEST_SIM_ROOT
source env.sh
cd $WEST_STRUCT_DATA_REF


# Make a temporary file in which to store output from the python script
#DIST=$(mktemp)


# Symlink a file needed for analysis
#ln -s $WEST_SIM_ROOT/namd_config/nacl.psf structure.psf


# Use a custom python script to calculate the distance between the Na+ and Cl-
# ions. This script looks for files named 'nacl.psf' and 'seg.dcd'.
source $WEST_SIM_ROOT/westpa_scripts/init_pcoord/jaccard.sh

while read i; do if [ "$i" = pcoord.txt ]; then break; fi; done

# Pipe the relevant part of the output file (the distance) to $WEST_PCOORD_RETURN
cp pcoord.txt $WEST_PCOORD_RETURN


if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
