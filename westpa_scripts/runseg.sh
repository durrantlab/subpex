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
ln -sv $WEST_SIM_ROOT/reference/mol.psf    structure.psf
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
# For multi-GPU on a single node, we need to explicitly assign the device consistent with
# the --n-workers=4.  4 is the max when the node has 4 GPUs.  Each worker will have its
# own thread id ($WM_PROCESS_INDEX) numbered from 0 to nGPU-1.

#export NAMD="$NAMD +devices $WM_PROCESS_INDEX"

echo "NAMD:"
echo $NAMD
echo "done"
$NAMD md.conf &> seg.log

########################## Calculate and return data ###########################

# Calculate the progress coordinate, which is the jaccard distance between the 
# bstate and this segment. 
# The script outputs the distance saving the values of the parent pcoord and the 
# child pcoord to a file called pcoord.txt.

cp -sv $WEST_SIM_ROOT/reference/mol.pdb ref.pdb
cp -sv $WEST_SIM_ROOT/westpa_scripts/settings.json .

###### Calculation of progress coordinate ######
#################### SubPEx ####################

# symlinks the parent pcoord.txt file and pipes the last line to the current pcoord
# file
ln -sv $WEST_PARENT_DATA_REF/pcoord.txt ./parentpcoord.txt
tail -n 1 parentpcoord.txt > pcoord.txt
ln -sv $WEST_PARENT_DATA_REF/pvol.txt ./parentpvol.txt
tail -n 1 parentpvol.txt > pvol.txt

# Check Chain for SubPEX Settings
#echo $(which python) >> /bgfs/jdurrant/durrantj/subpex_tests/NA/out
#echo $WEST_SIM_ROOT >> /bgfs/jdurrant/durrantj/subpex_tests/NA/out
#echo $(python $WEST_SIM_ROOT/westpa_scripts/jdistance.py ref.pdb seg.dcd settings.json) >> /bgfs/jdurrant/durrantj/subpex_tests/NA/out
#echo ====== >> /bgfs/jdurrant/durrantj/subpex_tests/NA/out

#python $WEST_SIM_ROOT/westpa_scripts/jdistance.py ref.pdb seg.dcd settings.json 2>&1 | grep -v "@>" >> pcoord.txt
python $WEST_SIM_ROOT/westpa_scripts/jdistance.py ref.pdb seg.dcd settings.json >> pcoord.txt

#python jdistance.py ../bstates/first_conformer/temp/ref.pdb ../reference/seg.dcd settings.json  2>&1 | grep -v "@>"

#paste <(cat jaccard.dat | awk {'print $2'}) <(cat rmsd.dat | awk {'print $2'}) > $WEST_PCOORD_RETURN



# this line just loops until we see the file 
while read i; do if [ -e timer.txt ]; then break; fi; done

#python2 $WEST_SIM_ROOT/westpa_scripts/test.py

cp pcoord.txt $WEST_PCOORD_RETURN
cp pvol.txt $WEST_PVOL_RETURN

# Clean up
rm -rf temp md.conf parent.coor parent.dcd parentpcoord.txt parentpvol.txt parent.vel parent.xsc ref.pdb settings.json timer.txt seg.coor.BAK seg.dcd.BAK seg.vel.BAK seg.xsc.BAK seg.xst.BAK
