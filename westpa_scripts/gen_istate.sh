#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

mkdir -p $(dirname $WEST_ISTATE_DATA_REF)
ln -s $WEST_BSTATE_DATA_REF $WEST_ISTATE_DATA_REF

mkdir $WEST_SIM_ROOT/reference
cp $WEST_SIM_ROOT/namd_config/mol.pdb $WEST_SIM_ROOT/reference/
cp $WEST_SIM_ROOT/namd_config/mol.psf $WEST_SIM_ROOT/reference/
cp $WEST_SIM_ROOT/prep/equil9/system-hc_eq_9.dcd  $WEST_SIM_ROOT/reference/seg.dcd