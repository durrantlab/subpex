from typing import Any

from collections.abc import Iterable

from loguru import logger

from .base import BaseContextManager, BaseContextValidator
from .cluster import ClusterContextManager
from .data import DataContextManager


class SubpexContextManager(
    BaseContextManager,
    DataContextManager,
    ClusterContextManager,
):
    """Contexts for SubPEx."""

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
        self.md_engine: str = "Amber"
        """Molecular dynamics engine to use. Both `Amber` and `NAMD` are implemented.
        """

        super().__init__(yaml_paths, **kwargs)


class SubpexContextValidator(BaseContextValidator):
    """Base class for validating Subpex contexts."""

    @staticmethod
    def write(value: Any, context: dict[str, Any]) -> bool:
        """Validate `write`"""
        return True

    @staticmethod
    def pocket_center(value: Any, context: dict[str, Any]) -> bool:
        """Validate `pocket_center`"""
        if len(value) != 3:
            logger.error("pocket_center must be a list of x, y, and z.")
            logger.error(f"pocket_center length should be 3, but is {len(value)}")
            return False
        for ele in value:
            if not isinstance(ele, float):
                logger.error(
                    f"pocket_center values must be `float` type; found {type(ele)}"
                )
                return False
        return True

    @staticmethod
    def pocket_radius(value: Any, context: dict[str, Any]) -> bool:
        """Validate `pocket_radius`"""
        if not isinstance(value, float):
            logger.error(f"pocket_radius value must be `float`; found {type(value)}")
            return False
        return True
