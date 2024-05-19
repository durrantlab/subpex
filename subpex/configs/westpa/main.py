from pydantic import BaseModel, Field

from ..io import YamlIO
from .data import DataConfig
from .executable import ExecutableConfig
from .propagation import PropagationConfig
from .system import SystemConfig


class WrapWestpaConfig(BaseModel):
    system: SystemConfig = Field(default_factory=SystemConfig)
    propagation: PropagationConfig = Field(default_factory=PropagationConfig)
    data: DataConfig = Field(default_factory=DataConfig)
    executable: ExecutableConfig = Field(default_factory=ExecutableConfig)


class WestpaConfig(BaseModel, YamlIO):
    """The root configuration for `west.cfg` files."""

    west: WrapWestpaConfig = Field(default_factory=WrapWestpaConfig)
