# Load the appropriate modules for mpi amber
#module purge

module load cuda/10.1
module load gcc/8.2.0 openmpi/4.0.3
module load amber/20
# Set up the Amber environment variable for mpi pmemd
export NPROC=$(nproc)
SANDER=pmemd.MPI
export AMBER="mpirun -n $NPROC $SANDER"

# Set the work manager to processes. Strange because NAMD itself is
# multithreaded, but it runs faster when each multitheraded NAMD
# is run simultaneously in its own thread.
export WORKMANAGER="processes"
