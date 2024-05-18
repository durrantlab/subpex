from typing import Any

from collections.abc import Sequence

import yaml
from pydantic import BaseModel


def write_yaml(config: BaseModel, file_name: str) -> None:
    """Serialize a Pydantic BaseModel instance to a YAML file.

    Args:
        config: The Pydantic BaseModel instance to be serialized.
        file_name: The name of the YAML file to write the serialized data to.

    Returns:
        None

    Raises:
        IOError: If the file cannot be written to.

    Example:

        ```python
        class MyConfig(BaseModel):
            param1: str
            param2: int

        config = MyConfig(param1="value1", param2=42)
        write_yaml(config, "config.yaml")
        ```
    """
    config_dict = config.dict()
    with open(file_name, "w") as f:
        yaml.dump(config_dict, f, default_flow_style=False)


def update_parameters_from_yaml(
    model_instance: Any, yaml_paths: Sequence[str] | str
) -> Any:
    """
    Update the attributes of a Pydantic BaseModel instance from one or more YAML files.

    Args:
        model_instance (Any): The Pydantic BaseModel instance to be updated.
        yaml_paths (Sequence[str] | str): A sequence of YAML file paths or a single YAML file path.

    Returns:
        Any: A new Pydantic BaseModel instance with updated attributes.

    Raises:
        FileNotFoundError: If any of the YAML files cannot be found.

    Example:
        ```python
        class MyConfig(BaseModel):
            param1: str
            param2: int

        config = MyConfig(param1="value1", param2=42)
        updated_config = update_parameters_from_yaml(config, "new_config.yaml")
        ```
    """
    if isinstance(yaml_paths, str):
        yaml_paths = [yaml_paths]
    for yaml_path in yaml_paths:
        with open(yaml_path, "r") as file:
            yaml_data = yaml.safe_load(file)
            model_instance = model_instance.copy(update=yaml_data)
    return model_instance
