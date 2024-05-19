from pydantic import Field
from pydantic_settings import BaseSettings

from .cluster import ClusteringConfig
from .data import DataConfig
from .pocket import PocketConfig


class SubpexConfig(BaseSettings):
    restart: bool = Field(default=False)
    """Restart a subpex simulation if directly already exists."""
    md_engine: str = Field(default="amber")
    """Molecular dynamics engine to use. Both `amber` and `namd` are supported."""

    pocket: PocketConfig = Field(default_factory=PocketConfig)
    data: DataConfig = Field(default_factory=DataConfig)
    clustering: ClusteringConfig = Field(default_factory=ClusteringConfig)
