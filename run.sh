#!/bin/bash
#
# run.sh
#
# Run the weighted ensemble simulation. Make sure to run init.sh first!
#

source env.sh

rm -f west.log

$WEST_ROOT/bin/w_run --work-manager $WORKMANAGER "$@" &> west.log
