#!/bin/bash
#
# get_pcoord.sh
#
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

# Use a custom script to calculate the jaccard distance between the starting 
# structure and the initial state (should be 0 since we are copying the files).
source $WEST_SIM_ROOT/westpa_scripts/pcoord_calc/jaccard.sh

# this line just loops until we see the file 
while read i; do if [ -e pcoord.txt ]; then break; fi; done

# Copy the file that contains the information to $WEST_PCOORD_RETURN. Better 
#than piping the info to the variable. 
cp pcoord.txt $WEST_PCOORD_RETURN

# Copy the file containing the volume of the pocket as calculated by subpex
cp pvol.txt $WEST_PVOL_RETURN

# If we are running in debug mode, then output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
