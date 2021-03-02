#!/bin/bash
#
# env.sh
#
# This script defines environment variables that are used by other shell
# scripts, both when setting up the simulation and when running the simulation.
#

module purge 
#KFW export PATH='/ihome/jdurrant/erh91/miniconda3/bin:$PATH'
export PATH="/ihome/jdurrant/erh91/miniconda3/bin:$PATH"
#. ~/.bashrc
#module load python/anaconda3.6-5.2.0
conda activate /ihome/jdurrant/erh91/miniconda3/envs/westpa-2020.02

################################## NAMD ########################################

# Set an environment variables for namd2. This is good practice for running
# WESTPA simulations on clusters, as repeatedly calling "namd2" or amber's pmemd rather than the 
# absolute path to the binaries (e.g., "/usr/local/namd/namd2") can take a long 
# time and be harder on the filesystem.  For this tutorial, this step is not 
# necessary but nonetheless demonstrates good practice.


source env-crc/mpi_namd.sh
#source env-crc/mpi_amber.sh

export NODELOC=$LOCAL
export USE_LOCAL_SCRATCH=1

############################## Python and WESTPA ###############################
# Next inform WESTPA what python it should use.  
#export WEST_PYTHON=$(which python2.7)
#source /ihome/jdurrant/erh91/westpa/westpa.sh
#echo $WEST_ROOT

# Check to make sure that the environment variable WEST_ROOT is set. 
# Here, the code '[ -z "$WEST_ROOT"]' will return TRUE if WEST_ROOT is not set,
# causing us to enter the if-statement, print an error message, and exit.
if [ -z "$WEST_ROOT" ]; then
  echo "The environment variable WEST_ROOT is not set."
  echo "Try running 'source westpa.sh' from the WESTPA installation directory"
  echo "This is going to cause problems unless you fix it!!!"
  #exit 1  # Because otherwise when you use source env.sh it leaves CRC!
fi

# Explicitly name our simulation root directory.  Similar to the statement 
# above, we check if the variable is not set.  If the variable is not set,
# the we set it 
#if [ -z "$WEST_SIM_ROOT" ]; then
#  export WEST_SIM_ROOT="$PWD"
#fi

export WEST_SIM_ROOT="$PWD"

# Set the simulation name.  Whereas "WEST_SIM_ROOT" gives us the entire 
# absolute path to the simulation directory, running the "basename" command
# will give us only the last part of that path (the directory name).
export SIM_NAME=$(basename $WEST_SIM_ROOT)

# WESTPA presumably uses it's own python installation. For processing
# trajectories, I'd like to use my python, which might have extra libraries
# like MDAnalysis installed.
#if [ -z "$MY_PYTHON" ]; then
#  echo "The environment variable MY_PYTHON is not set."
#  echo "I'll just use the same python WESTPA is using."
#  export MY_PYTHON=$(which python2.7)
#fi

export WM_ZMQ_MASTER_HEARTBEAT=100
export WM_ZMQ_WORKER_HEARTBEAT=100
export WM_ZMQ_TIMEOUT_FACTOR=100
export BASH=$SWROOT/bin/bash
export PERL=$SWROOT/usr/bin/perl
export ZSH=$SWROOT/bin/zsh
export IFCONFIG=$SWROOT/bin/ifconfig
export CUT=$SWROOT/usr/bin/cut
export TR=$SWROOT/usr/bin/tr
export LN=$SWROOT/bin/ln
export CP=$SWROOT/bin/cp
export RM=$SWROOT/bin/rm
export SED=$SWROOT/bin/sed
export CAT=$SWROOT/bin/cat
export HEAD=$SWROOT/bin/head
export TAR=$SWROOT/bin/tar
export AWK=$SWROOT/usr/bin/awk
export PASTE=$SWROOT/usr/bin/paste
export GREP=$SWROOT/bin/grep
export SORT=$SWROOT/usr/bin/sort
export UNIQ=$SWROOT/usr/bin/uniq
export HEAD=$SWROOT/usr/bin/head
export MKDIR=$SWROOT/bin/mkdir
export ECHO=$SWROOT/bin/echo
export DATE=$SWROOT/bin/date