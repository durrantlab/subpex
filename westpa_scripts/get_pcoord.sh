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

rm -rf temp
mkdir temp

cp $WEST_STRUCT_DATA_REF/seg.dcd    temp/
cp $WEST_SIM_ROOT/reference/mol.pdb temp/ref.pdb
cp $WEST_SIM_ROOT/westpa_scripts/settings.json temp/

# Use a custom script to calculate the jaccard distance between the starting 
# structure and the initial state (should be 0 since we are copying the files).
python $WEST_SIM_ROOT/westpa_scripts/jdistance.py temp/ref.pdb temp/seg.dcd temp/settings.json > pcoord.txt

#paste <(cat jaccard.dat | awk {'print $2'}) <(cat rmsd.dat | awk {'print $2'}) > $WEST_PCOORD_RETURN

# this line just loops until we see the file 
while read i; do if [ -e pcoord.txt ]; then break; fi; done

cp pcoord.txt $WEST_PCOORD_RETURN

# Copy the file that contains the information to $WEST_PCOORD_RETURN. Better 
#than piping the info to the variable. 
#cp pcoord.txt $WEST_PCOORD_RETURN
#tail -n 1 pcoord.txt > $WEST_PCOORD_RETURN

# Copy the file containing the volume of the pocket as calculated by subpex
tail -n 1 pvol.txt > $WEST_PVOL_RETURN

# If we are running in debug mode, then output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
