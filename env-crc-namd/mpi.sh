# Load the appropriate modules for mpi namd
#module purge

if [ $MD_ENGINE == "NAMD" ]
then
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
elif [ $MD_ENGINE == "AMBER" ]
then 
    module load cuda/10.1
    module load gcc/8.2.0 openmpi/4.0.3
    module load amber/20
    # Set up the NAMD environment variable for mpi namd
    export NPROC=$(nproc)
    SANDER=pmemd.MPI
    export AMBER="mpirun -n $NPROC $SANDER"

    # Set the work manager to processes. Strange because NAMD itself is
    # multithreaded, but it runs faster when each multitheraded NAMD
    # is run simultaneously in its own thread.
    export WORKMANAGER="processes"
else
    echo "This engine is not suported"
fi