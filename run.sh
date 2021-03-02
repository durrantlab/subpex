#!/bin/bash
#
# run.sh
#
# Run the weighted ensemble simulation. Make sure you ran init.sh first!
#

source env.sh

rm -f west.log

# For multiple threads. But not accross nodes. For more, see:
# https://github.com/westpa/westpa/wiki/Running-WESTPA-in-a-multi-node-environment
$WEST_ROOT/bin/w_run --work-manager $WORKMANAGER "$@" &> west.log


