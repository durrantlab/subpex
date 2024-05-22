from pydantic import BaseModel, Field

from ..io import YamlIO
from .data import DataConfig
from .environment import WestpaEnv
from .executable import ExecutableConfig
from .propagation import PropagationConfig
from .system import SystemConfig


class WrapWestpaConfig(BaseModel):
    """
    Configuration class for wrapping WESTPA configuration options.

    Attributes:
        system (SystemConfig): Configuration options for the system.
        propagation (PropagationConfig): Configuration options for the propagation.
        data (DataConfig): Configuration options for the data.
        executable (ExecutableConfig): Configuration options for the executable.
    """

    system: SystemConfig = Field(default_factory=SystemConfig)
    propagation: PropagationConfig = Field(default_factory=PropagationConfig)
    data: DataConfig = Field(default_factory=DataConfig)
    executable: ExecutableConfig = Field(default_factory=ExecutableConfig)


class WestpaConfig(BaseModel, YamlIO):
    """The root configuration for `west.cfg` files."""

    west: WrapWestpaConfig = Field(default_factory=WrapWestpaConfig)
    env: WestpaEnv = Field(default_factory=WestpaEnv)
