from collections.abc import MutableMapping, Sequence

from pydantic import BaseModel


class SlurmConfig(BaseModel):
    """Context manager for Slurm job submission scripts.

    This class provides a structured way to define and manage the configuration for
    submitting jobs to a Slurm workload manager. Each attribute corresponds to a
    specific Slurm configuration parameter or job setup step.
    """

    job_name: str = "subpex_job"
    """Unique name for the Slurm job.

    This should be descriptive to help identify the job among many others in the queue.
    Example: "data_analysis_job"
    """

    nodes: int = 1
    """Number of nodes to use for the Slurm job.

    Adjust this based on the job's resource requirements. For instance, a large parallel
    job might need several nodes, while a smaller job might only need one.
    Example: 4
    """

    ntasks_per_node: int = 8
    """Number of tasks to run per node.

    This typically corresponds to the number of CPU cores to use on each node. Adjust this
    based on the node's capabilities and the parallelism of your job.
    Example: 16
    """

    output_path: str = "slurm-%j.out"
    """Path for the job's standard output file.

    Use `%j` to include the job ID in the filename, ensuring that each job's output is
    saved to a unique file.
    Example: "logs/job_output_%j.out"
    """

    error_path: str = "slurm-%j.err"
    """Path for the job's error output file.

    Similar to `output_path`, using `%j` in the filename ensures that errors for each job
    are logged separately.
    Example: "logs/job_errors_%j.err"
    """

    time: str = "1-00:00:00"
    """Maximum time for the job.

    Specified in the format `D-HH:MM:SS`. Adjust this based on the expected runtime of your job.
    Example: "0-12:00:00" for a 12-hour job.
    """

    cluster: str = "smp"
    """Cluster name where the job will run.

    Ensure this matches the available cluster names in your Slurm environment. This helps
    direct the job to the appropriate set of resources.
    Example: "hpc_cluster"
    """

    partition: str = "smp"
    """Partition name to submit the job to.

    Choose an appropriate partition based on resource needs and availability. Partitions
    can have different resource limits and policies.
    Example: "short"
    """

    account: str | None = None
    """Account name for job submission.

    Set this to your project's account name to properly attribute resource usage, or leave
    it as None if not needed.
    Example: "research_project_123"
    """

    modules: Sequence[str] = []
    """List of modules to load before running the job.

    Include all necessary software modules that your job requires. This ensures the
    environment is correctly set up before execution.
    Example: ["python/3.8", "gcc/9.2"]
    """

    env_vars: MutableMapping[str, str] = {}
    """Dictionary of environment variables to set before running the job.

    Use this to configure the job's environment, setting any necessary environment
    variables.
    Example: {"OMP_NUM_THREADS": "16", "MY_VARIABLE": "value"}
    """

    commands_pre: Sequence[str] = []
    """List of commands to run before the main job command.

    Useful for setup tasks like copying files, creating directories, or loading
    additional software. These commands will be executed before the main job starts.
    Example: ["mkdir -p /scratch/my_job", "cp input.dat /scratch/my_job/"]
    """

    commands_run: Sequence[str] = []
    """List of main commands to run for the job.

    This should include the primary executable or script for the job. These are the
    main tasks that the job will perform.
    Example: ["python my_script.py", "./run_simulation.sh"]
    """

    commands_post: Sequence[str] = []
    """List of commands to run after the main job command.

    Use this for cleanup tasks or additional processing. These commands will be executed
    after the main job tasks are completed.
    Example: ["cp /scratch/my_job/output.dat .", "rm -rf /scratch/my_job"]
    """
