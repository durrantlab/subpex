from typing import Any

from collections.abc import Iterable, MutableMapping

from loguru import logger
from ruamel.yaml import YAML

from .cluster import ClusterContextManager
from .data import DataContextManager
from .pocket import PocketContextManager
from .slurm import SlurmContextManager


class SubpexContextManager(
    DataContextManager,
    PocketContextManager,
    ClusterContextManager,
    SlurmContextManager,
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

        DataContextManager.__init__(self, yaml_paths, **kwargs)
        PocketContextManager.__init__(self, yaml_paths, **kwargs)
        ClusterContextManager.__init__(self, yaml_paths, **kwargs)
        SlurmContextManager.__init__(self, yaml_paths, **kwargs)

        self.yaml_path: str | None = None
        if isinstance(yaml_paths, str):
            yaml_paths = [yaml_paths]
        if yaml_paths is not None:
            for yaml_path in yaml_paths:
                self.from_yaml(yaml_path)
        self.update(kwargs)

    def from_yaml(self, yaml_path: str | None) -> None:
        """Load context information from a YAML file. This will only update data
        contained in the YAML file.

        Args:
            yaml_path: Path to YAML file to load.
        """
        if yaml_path is not None:
            logger.info("Loading YAML context from {}", yaml_path)
            yaml = YAML(typ="safe")
            with open(yaml_path, "r", encoding="utf-8") as f:
                yaml_data = yaml.load(f)
            logger.debug("YAML data:\n{}", yaml_data)
            self.update(yaml_data)
        self.yaml_path = yaml_path

    def update(self, attr_dict: MutableMapping[str, Any]) -> dict[str, Any]:
        """Update attributes with values from the provided dictionary.

        Args:
            attr_dict: Dictionary containing attribute names and their
            corresponding values.
        """
        logger.debug("Updating context:\n{}", attr_dict)
        for key, value in attr_dict.items():
            setattr(self, key, value)
        return self.get()

    def get(self) -> dict[str, Any]:
        """Retrieve the context.

        Returns:
            A dictionary representing the current context.
        """
        # The following line filters methods and attributes like __dict__.
        context = {
            k: v for k, v in vars(self).items() if not callable(v) and "__" not in k
        }
        logger.debug("Retrieved context:\n{}", context)
        return context

    def __enter__(self) -> dict[str, Any]:
        """Enter the context and return the current context as a dictionary."""
        return self.get()

    def __exit__(self, exc_type, exc_value, exc_tb):
        """Exit the context.

        Args:
            exc_type: Type of the exception.
            exc_value: Value of the exception.
            exc_tb: Traceback information.
        """


class SubpexContextValidator:
    """Base class for validating Subpex contexts."""

    @classmethod
    def validate(cls, context_manager: SubpexContextManager) -> bool:
        """Validate contexts for subpex.

        Args:
            context_manager: A subpex context manager to validate.

        Returns:
            If the context is valid.
        """
        logger.info("Validating with: {}", cls.__name__)
        is_valid = True
        context = context_manager.get()
        for key, value in context.items():
            try:
                checker = getattr(cls, key)
            except AttributeError:
                logger.debug("Cannot check {}", key)
                continue
            logger.debug("Checking {}", key)
            if value is not None:
                if isinstance(checker, tuple):
                    if value not in checker:
                        logger.error("  Invalid: {}", value)
                        is_valid = False
                    else:
                        logger.debug("  Valid: {}", value)
                if callable(checker):
                    is_value_valid = checker(value, context)
                    if not is_value_valid:
                        logger.error("  Invalid: {}", value)
                        is_valid = False
                    else:
                        logger.debug("  Valid: {}", value)
            else:
                logger.debug("  Skipping: None")
        return is_valid

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
