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
#     directory for running pmemd/sander oe run NAMD2
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

# Either continue an existing tractory, or start a new trajectory. Here, both
# cases are the same.  If you need to handle the cases separately, you can
# check the value of the environment variable "WEST_CURRENT_SEG_INIT_POINT",
# which is equal to either "SEG_INITPOINT_CONTINUES" or "SEG_INITPOINT_NEWTRAJ"
# for continuations of previous segments and new trajectories, respecitvely.
# For an example, see the nacl_amb tutorial.

# The weighted ensemble algorithm requires that dynamics are stochastic.
# We'll use the "sed" command to replace the string "RAND" with a randomly
# generated seed.

# This trajectory segment will start off where its parent segment left off.
# The "ln" command makes symbolic links to the parent segment's edr, gro, and 
# and trr files. This is preferable to copying the files, since it doesn't
# require writing all the data again.
ln -sv $WEST_PARENT_DATA_REF/pcoord.txt ./parent_pcoord.txt
[[ -e $WEST_PARENT_DATA_REF/rog.txt ]] && ln -sv $WEST_PARENT_DATA_REF/rog.txt ./parent_rog.txt
[[ -e $WEST_PARENT_DATA_REF/bb.txt ]] && ln -sv $WEST_PARENT_DATA_REF/bb.txt ./parent_bb.txt
[[ -e $WEST_PARENT_DATA_REF/pvol.txt ]] && ln -sv $WEST_PARENT_DATA_REF/pvol.txt ./parent_pvol.txt
[[ -e $WEST_PARENT_DATA_REF/jd.txt ]] && ln -sv $WEST_PARENT_DATA_REF/jd.txt ./parent_jd.txt
[[ -e $WEST_PARENT_DATA_REF/prmsd.txt ]] && ln -sv $WEST_PARENT_DATA_REF/prmsd.txt ./parent_prmsd.txt
[[ -e $WEST_PARENT_DATA_REF/fop.txt ]] && ln -sv $WEST_PARENT_DATA_REF/fop.txt ./parent_fop.txt
[[ -e $WEST_PARENT_DATA_REF/composite.txt ]] && ln -sv $WEST_PARENT_DATA_REF/composite.txt ./parent_composite.txt

# CODE TO SET RANDOM SEED

# IF USING THE NAMD MD ENGINE:
sed "s/RAND/$WEST_RAND16/g" \
$WEST_SIM_ROOT/reference/namd.md.conf > namd.md.conf

# # IF USING THE AMBER MD ENGINE:
# sed "s/RAND/$WEST_RAND16/g" \
# $WEST_SIM_ROOT/reference/amber.prod_npt.in > amber.prod_npt.in

# This trajectory segment will start off where its parent segment left off.
# The "ln" command makes symbolic links to the parent segment's edr, gro, and 
# and trr files. This is preferable to copying the files, since it doesn't
# require writing all the data again.

# LINK FILES TO RESTART WALKER

# NAMD files
# ln -sv $WEST_SIM_ROOT/reference/mol.prmtop mol.prmtop
# ln -sv $WEST_SIM_ROOT/reference/mol.inpcrd mol.inpcrd
# ln -sv $WEST_PARENT_DATA_REF/seg.coor ./parent.coor
# ln -sv $WEST_PARENT_DATA_REF/seg.dcd  ./parent.dcd
# ln -sv $WEST_PARENT_DATA_REF/seg.vel  ./parent.vel
# ln -sv $WEST_PARENT_DATA_REF/seg.xsc  ./parent.xsc

# AMBER files
ln -sv $WEST_PARENT_DATA_REF/seg.rst ./parent.rst

############################## Run the dynamics ################################

# CALL MD ENGINE

# NAMD MD engine
$NAMD namd.md.conf > seg.log

if grep -q RATTLE seg.log; then
    $NAMD namd.md.conf > seg.log
fi

# # Amber MD engine
# $AMBER -O -i amber.prod_npt.in -p mol.prmtop -c parent.rst -r seg.rst -x seg.nc -o seg.log -inf seg.nfo

########################## Calculate and return progress coordiante ###########################
######################################### SubPEx ##############################################

# Calculate the progress coordinate, which is the jaccard distance between the 
# bstate and this segment. 
# The script outputs the distance saving the values of the parent pcoord and the 
# child pcoord to a file called pcoord.txt.

# RUN PCOORD.PY

# FOR NAMD
python3 $WEST_SIM_ROOT/westpa_scripts/pcoord.py seg.dcd  $WEST_SIM_ROOT/west.cfg --we

# # FOR AMBER:
# python3 $WEST_SIM_ROOT/westpa_scripts/pcoord.py seg.nc  $WEST_SIM_ROOT/west.cfg --we

cp pcoord.txt $WEST_PCOORD_RETURN
[[ -e pvol.txt ]] && cp pvol.txt $WEST_PVOL_RETURN
[[ -e rog.txt ]] && cp rog.txt $WEST_ROG_RETURN
[[ -e bb.txt ]] && cp bb.txt $WEST_BB_RETURN
[[ -e jd.txt ]] && cp jd.txt $WEST_JD_RETURN
[[ -e prmsd.txt ]] && cp prmsd.txt $WEST_PRMSD_RETURN
[[ -e fop.txt ]] && cp fop.txt $WEST_FOP_RETURN  # TODO: Erich, check this.
[[ -e composite.txt ]] && cp composite.txt $WEST_COMPOSITE_RETURN

# Clean up
if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
else
  rm parent.* mol.prmtop
  #rm timer.txt seg.coor.BAK seg.dcd.BAK seg.vel.BAK seg.xsc.BAK seg.xst.BAK
  #rm *restart*
fi

