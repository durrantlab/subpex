#!/bin/bash

#### Prepare Temporary Space ####
rm -rf temp
mkdir temp

SCRIPTS=$WEST_SIM_ROOT/westpa_scripts/pcoord_calc/

#################### Trajectory to PDB ####################

cp $WEST_STRUCT_DATA_REF/seg.dcd    temp/
cp $WEST_SIM_ROOT/reference/mol.psf temp/
cp $WEST_SIM_ROOT/reference/mol.pdb temp/ref.pdb
cp $SCRIPTS/SubPEX_settings_pre     ./

######
python2 $SCRIPTS/align_traj.py temp/mol.psf seg.dcd temp/ref.pdb

#################### SubPEX ####################

# Check Chain for SubPEX Settings
python2 $SCRIPTS/fix_settings.py temp/seg_aligned.pdb

# Run SubPEX
python2 $SCRIPTS/SubPEX_tweaked.py


#################### Jaccard Index #################### 

# Make & Copy Reference
cp temp/seg_aligned_frame0.xyz $WEST_STRUCT_DATA_REF/ref.xyz
cp temp/seg_aligned_frame0.xyz $WEST_SIM_ROOT/reference/ref.xyz

# Run Conversion
python2 $SCRIPTS/jdistance.py > pcoord.txt

#### Delete Temporary Space ####
#rm -rf temp
