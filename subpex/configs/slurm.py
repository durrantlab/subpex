from typing import Any

from collections.abc import Iterable, MutableMapping, Sequence


class SlurmConfig:
    """Context manager for Slurm job submission scripts."""

    def __init__(
        self, yaml_paths: str | Iterable[str] | None = None, **kwargs: dict[str, Any]
    ) -> None:
        """
        Args:
            yaml_paths: Path(s) to YAML file(s) to load into the context.
            kwargs: Additional keyword arguments to initialize the context.
        """
        self.job_name: str = "subpex_job"
        """Unique name for the Slurm job. This should be descriptive to help identify
        the job."""

        self.nodes: int = 1
        """Number of nodes to use for the Slurm job. Adjust this based on the job's
        resource requirements.
        """

        self.ntasks_per_node: int = 8
        """Number of tasks to run per node. This typically corresponds to the number of
        CPU cores to use on each node.
        """

        self.output_path: str = "slurm-%j.out"
        """Path for the job's standard output file. Use %j to include the job ID in the
        filename."""

        self.error_path: str = "slurm-%j.err"
        """Path for the job's error output file. Use %j to include the job ID in the
        filename."""

        self.time: str = "1-00:00:00"
        """Maximum time for the job in the format D-HH:MM:SS. Adjust based on expected runtime."""

        self.cluster: str = "smp"
        """Cluster name where the job will run. Ensure this matches the available cluster names."""

        self.partition: str = "smp"
        """Partition name to submit the job to. Choose an appropriate partition based on resource needs and availability."""

        self.account: str | None = None
        """Account name for job submission. Set this to your project's account name, or leave it as None if not needed."""

        self.modules: Sequence[str] = []
        """List of modules to load before running the job. Include all necessary software modules."""

        self.env_vars: MutableMapping[str, str] = {}
        """Dictionary of environment variables to set before running the job. Use this to configure the job's environment."""

        self.commands_pre: Sequence[str] = []
        """List of commands to run before the main job command. Useful for setup tasks like copying files or creating directories."""

        self.commands_run: Sequence[str] = []
        """List of main commands to run for the job. This should include the primary executable or script for the job."""

        self.commands_post: Sequence[str] = []
        """List of commands to run after the main job command. Use this for cleanup tasks or additional processing."""
