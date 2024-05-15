from typing import Any

from collections.abc import Iterable

from loguru import logger

from .base import BaseContextManager, BaseContextValidator


class PCoordContextManager(BaseContextManager):
    """Contexts for progress coordinate."""

    # pylint: disable-next=too-many-statements
    def __init__(
        self, yaml_paths: str | Iterable[str] | None = None, **kwargs: dict[str, Any]
    ) -> None:
        """
        Args:
            yaml_paths: Path(s) to YAML file(s) to load into the context.
        """
        self.restart: bool = False
        """Restart a subpex simulation if directly already exists.

        TODO: run the sp_restart script if westpa file exists and this is True.
        """

        super().__init__(yaml_paths, **kwargs)


# pylint: disable-next=too-few-public-methods
class PCoordContextValidator(BaseContextValidator):
    """Base class for validating contexts."""

    # pylint: disable=unused-argument
