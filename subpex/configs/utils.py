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
