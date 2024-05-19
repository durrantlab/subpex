from pydantic import BaseModel, Field

from .system import System


class WrapWestpaConfig(BaseModel):
    system: System = Field(default_factory=System)


class WestpaConfig(BaseModel):
    """The root configuration for `west.cfg` files."""

    west: WrapWestpaConfig = Field(default_factory=WrapWestpaConfig)
