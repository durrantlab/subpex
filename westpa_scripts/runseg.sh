#!/bin/bash
#
# runseg.sh
#
# WESTPA runs this script for each trajectory segment. WESTPA supplies
# environment variables that are unique to each segment, such as:
#
#   WEST_CURRENT_SEG_DATA_REF: A path to where the current trajectory segment's
#       data will be stored. This will become "WEST_PARENT_DATA_REF" for any
#       child segments that spawn from this segment
#   WEST_PARENT_DATA_REF: A path to a file or directory containing data for the
#       parent segment.
#   WEST_CURRENT_SEG_INITPOINT_TYPE: Specifies whether this segment is starting
#       anew, or if this segment continues from where another segment left off.
#   WEST_RAND16: A random integer
#
# This script has the following three jobs:
#  1. Create a directory for the current trajectory segment, and set up the
#     directory for running pmemd/sander 
#  2. Run the dynamics
#  3. Calculate the progress coordinates and return data to WESTPA


# If we are running in debug mode, then output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

######################## Set up for running the dynamics #######################

# To be sure Load environment
source $WEST_SIM_ROOT/env.sh

# Set up the directory where data for this segment will be stored.
cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

# Make symbolic links to the topology file and parameter files. These are not 
# unique to each segment.
ln -sv $WEST_SIM_ROOT/reference/mol.pdb    structure.pdb
ln -sv $WEST_SIM_ROOT/reference/mol.prmtop structure.prmtop
ln -sv $WEST_SIM_ROOT/reference/mol.inpcrd structure.inpcrd

# Either continue an existing tractory, or start a new trajectory. Here, both
# cases are the same.  If you need to handle the cases separately, you can
# check the value of the environment variable "WEST_CURRENT_SEG_INIT_POINT",
# which is equal to either "SEG_INITPOINT_CONTINUES" or "SEG_INITPOINT_NEWTRAJ"
# for continuations of previous segments and new trajectories, respecitvely.
# For an example, see the nacl_amb tutorial.

# The weighted ensemble algorithm requires that dynamics are stochastic.
# We'll use the "sed" command to replace the string "RAND" with a randomly
# generated seed.
sed "s/RAND/$WEST_RAND16/g" \
  $WEST_SIM_ROOT/reference/md.conf > md.conf

# This trajectory segment will start off where its parent segment left off.
# The "ln" command makes symbolic links to the parent segment's edr, gro, and 
# and trr files. This is preferable to copying the files, since it doesn't
# require writing all the data again.
ln -sv $WEST_PARENT_DATA_REF/seg.coor ./parent.coor
ln -sv $WEST_PARENT_DATA_REF/seg.dcd  ./parent.dcd
ln -sv $WEST_PARENT_DATA_REF/seg.vel  ./parent.vel
ln -sv $WEST_PARENT_DATA_REF/seg.xsc  ./parent.xsc

############################## Run the dynamics ################################
# Propagate the segment using namd2 
$NAMD md.conf > seg.log 

if grep -q RATTLE seg.log; then
    $NAMD md.conf > seg.log
fi

########################## Calculate and return data ###########################

# Calculate the progress coordinate, which is the jaccard distance between the 
# bstate and this segment. 
# The script outputs the distance saving the values of the parent pcoord and the 
# child pcoord to a file called pcoord.txt.

cp -sv $WEST_SIM_ROOT/reference/mol.pdb ref.pdb
cp -sv $WEST_SIM_ROOT/westpa_scripts/settings.json .

###### Calculation of progress coordinate ######
#################### SubPEx ####################

if [ -n "$SEG_DEBUG" ] ; then
  python3 $WEST_SIM_ROOT/westpa_scripts/jdistance.py ref.pdb structure.prmtop seg.dcd settings.json --debug
else
  python3 $WEST_SIM_ROOT/westpa_scripts/jdistance.py ref.pdb structure.prmtop seg.dcd settings.json
fi

cp pcoord.txt $WEST_PCOORD_RETURN
cp pvol.txt $WEST_PVOL_RETURN
cp rog.txt $WEST_ROG_RETURN
cp bb_rmsd.txt $WEST_BB_RETURN

# Clean up
if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
else
  rm parent.coor parent.dcd  parent.vel parent.xsc structure.inpcrd structure.prmtop
  rm ref.pdb settings.json timer.txt seg.coor.BAK seg.dcd.BAK seg.vel.BAK seg.xsc.BAK seg.xst.BAK
  rm pcoord.txt pvol.txt rog.txt bb_rmsd.txt seg.restart.coor seg.restart.coor.old seg.restart.vel 
  rm seg.restart.vel.old seg.restart.xsc seg.restart.xsc.old
fi

