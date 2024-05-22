from pydantic import BaseModel, Field


class WestpaEnv(BaseModel):
    """Configuration for WESTPA environment variables.

    This is used by
    [write_env_script][setup.westpa.scripts.write_env_script]
    to create the `env.sh` script for setting up the WESTPA environment.
    """

    WEST_SIM_ROOT: str = Field(default="$PWD")
    """path to the base directory containing the WESTPA install"""
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
