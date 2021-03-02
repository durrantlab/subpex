#!/bin/bash
#
# env.sh
#
# This script defines environment variables that are used by other shell
# scripts, both when setting up the simulation and when running the simulation.
#



# Check to make sure that the environment variable WEST_ROOT is set. 
# Here, the code '[ -z "$WEST_ROOT"]' will return TRUE if WEST_ROOT is not set,
# causing us to enter the if-statement, print an error message, and exit.
if [ -z "$WEST_ROOT" ]; then
  echo "The environment variable WEST_ROOT is not set."
  echo "Try running 'source westpa.sh' from the WESTPA installation directory"
  echo "This is going to cause problems unless you fix it!!!"
fi

# Explicitly name our simulation root directory.  Similar to the statement 
# above, we check if the variable is not set.  If the variable is not set,
# the we set it 
if [ -z "$WEST_SIM_ROOT" ]; then
  export WEST_SIM_ROOT="$PWD"
fi

export WEST_SIM_ROOT="$PWD"

# Set the simulation name.  Whereas "WEST_SIM_ROOT" gives us the entire 
# absolute path to the simulation directory, running the "basename" command
# will give us only the last part of that path (the directory name).
export SIM_NAME=$(basename $WEST_SIM_ROOT)




