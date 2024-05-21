from pydantic import Field
from pydantic_settings import BaseSettings

from ..io import YamlIO
from .cluster import ClusteringConfig
from .data import DataConfig
from .pocket import PocketConfig
from .sim import SimulationConfig


class SubpexConfig(BaseSettings, YamlIO):
    restart: bool = Field(default=False)
    """Restart a subpex simulation if directly already exists."""
    calculated_points: int = Field(default=-1)
    """Number of point to calculate per trajectory segment. If `-1`, it will
    calculate all.
    """

    pocket: PocketConfig = Field(default_factory=PocketConfig)
    data: DataConfig = Field(default_factory=DataConfig)
    clustering: ClusteringConfig = Field(default_factory=ClusteringConfig)
    sim: SimulationConfig = Field(default_factory=SimulationConfig)
