from collections.abc import MutableSequence

from pydantic import BaseModel, Field


class EnvVars(BaseModel):
    WEST_SIM_ROOT: str = Field(default="$PWD")
    """path to the base directory containing the WESTPA install"""
    SIM_NAME: str = Field(default="$(basename $WEST_SIM_ROOT)")
    """Simulation name"""
    WEST_PYTHON: str = Field(default="$(which python)")
    """path to python executable to run the WESTPA simulation"""
    WEST_PYTHONPATH: str | None = Field(default=None)
    """path to any additional modules that WESTPA will require to run the simulation"""
    WEST_CURRENT_ITER: int = Field(ge=1, default=1)
    """Current iteration number"""
    WEST_CURRENT_SEG_ID: int = Field(ge=0, default=0)
    """Current segment ID"""
    WEST_CURRENT_SEG_DATA_REF: str | None = Field(default=None)
    """General-purpose reference, based on current segment information, configured in west.cfg. Usually used for storage paths. """
    WEST_PARENT_ID: int | None = Field(default=None)
    """Segment ID of parent segment. Negative for initial points."""
    WEST_PCOORD_RETURN: str | None = Field(default=None)
    """Where progress coordinate data must be stored."""


class WestpaEnv(BaseModel):
    """Configuration for WESTPA environment setup script; typically called `env.sh`.

    This is used by
    [write_env_script][setup.westpa.scripts.write_env_script]
    to create the `env.sh` script for setting up the WESTPA environment.
    """

    vars: EnvVars = Field(default_factory=EnvVars)
    """Keys and values for the environment variables to be set in the `env.sh` script.

    Additional environment variables can be added with the following syntax:

    ```python
    WestpaConfig.env.vars.KEY = VALUE
    ```

    For example, setting the number of processors (`NPROCS`) to `12` would be

    ```python
    WestpaConfig.env.vars.NPROCS = 12
    ```
    """

    purge_modules: bool = Field(default=True)
    """If True, remove all loaded modules before loading the modules in `load_modules`."""

    load_modules: MutableSequence[str] = Field(default=[])
    """List of modules to load before running WESTPA simulation.
    Each module is loaded using `module load` command.

    Warning:
        These modules are specific to the cluster and the software version and often
        need to be loaded in a specific order. Please consult the cluster documentation.

    Example:
        If we were to use AMBER 22 for a WESTPA simulation on the
        University of Pittsburgh's H2P cluster, we would set this to

        ```python
        load_modules = ["gcc/10.2.0", "openmpi/4.1.1", "amber/22"]
        ```

        For NAMD 2.13, it would be

        ```python
        load_modules = ["intel/2018.2.199", "openmpi/2.1.1", "namd/3.1.1"]
        ```

        For GROMACS 2021.2 it would be

        ```python
        load_modules = ["gcc/10.2.0", "openmpi/4.1.1", "gromacs/2021.2"]
        ```
    """

    lines_prepend: MutableSequence[str] = Field(default=[])
    """Lines to prepend to the `env.sh` script."""
    lines_append: MutableSequence[str] = Field(default=[])
    """Lines to append to the `env.sh` script."""
