from pydantic import BaseModel, Field

from .data import DataConfig
from .executable import ExecutableConfig
from .propagation import PropagationConfig
from .system import SystemConfig


class WestpaConfigConfig(BaseModel):
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
