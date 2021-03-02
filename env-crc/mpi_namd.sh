# Load the appropriate modules for mpi namd
#module purge

module load intel/2017.1.132 intel-mpi/2017.1.132
module load namd/2.12

# Set up the NAMD environment variable for mpi namd
export NPROC=$(nproc)
export CHARMRUN=$(which charmrun)
export NAMDBIN=$(which namd2)
export NAMD="$CHARMRUN +p$NPROC $NAMDBIN"

# Set the work manager to processes. Strange because NAMD itself is
# multithreaded, but it runs faster when each multitheraded NAMD
# is run simultaneously in its own thread.
export WORKMANAGER="processes"
