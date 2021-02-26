# Load the appropriate modules for mpi namd
module purge
module load namd/2.12b1-multicore-CUDA

# Set up the NAMD environment variable for mpi namd
#KFW export NPROC=$(nproc)
export NPROC=12
export NAMDBIN=$(which namd2)
export NAMD="$NAMDBIN +p$NPROC +idlepoll +setcpuaffinity"

#KFW #export NAMD="$NAMDBIN +p$NPROC +idlepoll +setcpuaffinity"
#KFW export NAMD="$NAMDBIN +p12 +idlepoll +setcpuaffinity"

# Useful note from Kim re. above line:
# Yes, strictly speaking the threads should be physical CPU cores/ GPUs but sometime you can over subscript the number of physical CPU cores and still get some benefit.  You have to use trial and error to determine how much you can exceed the number of physical cores.

# Set the work manager to processes. Strange because NAMD itself is
# multithreaded, but it runs faster when each multitheraded NAMD
# is run simultaneously in its own thread.
#export WORKMANAGER="processes"
export WORKMANAGER="serial"
