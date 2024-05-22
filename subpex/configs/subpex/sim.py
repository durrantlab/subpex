from pydantic import BaseModel, Field


class SimulationConfig(BaseModel):

    md_engine: str = Field(default="amber")
    """Molecular dynamics engine to use. Both `amber` and `namd` are supported."""
    path_traj_relax: str | None = Field(default=None)
    """Path to the final simulation trajectory file to use. The last frame of this
    trajectory will be used as the starting frame for the enhanced simulations.
    """
    path_topo: str | None = Field(default=None)
    """Path to the topology file for the system. This file is used to define the system
    and parameters for the molecular dynamics simulations.
    """
