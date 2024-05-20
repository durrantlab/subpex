import argparse
import os
from collections.abc import Sequence

from loguru import logger

from .main import WestpaConfig


def run_westpa_config(
    file_name: str = "west.cfg",
    save_dir: str | None = None,
    yaml_paths: Sequence[str] | None = None,
) -> None:
    """
    Create, update, and save a WESTPA configuration file.

    Args:
        file_name: The name of the output configuration file.
        save_dir: The directory where the configuration file will be saved. If `None`,
            the current directory is used.
        yaml_paths: A list of paths to YAML files from which to update the
            configuration. If None, no updates are made.

    Returns:
    --------
    None
    """
    logger.debug("Initializing default WestpaConfig instance.")
    westpa_config = WestpaConfig()

    if yaml_paths:
        logger.debug("Updating parameters from YAML files.")
        westpa_config.from_yaml(yaml_paths)

    output_path = os.path.join(save_dir or "", file_name)
    logger.debug(f"Writing the YAML configuration to {output_path}.")
    westpa_config.to_yaml(output_path)
    logger.debug("Configuration successfully saved.")


def cli_westpa_config():
    """Command-line interface for generating a WESTPA configuration file.

    This function provides a CLI for the `run_westpa_config` function, allowing users
    to generate and save a WESTPA configuration file by specifying the file name,
    save directory, and paths to YAML files for updating the configuration parameters.

    Args:
        --name: The name of the output configuration file including the file extension.
        --save_dir: The directory where the configuration file will be saved.
        --yaml: One or more paths to YAML files from which to update the configuration
            parameters. These are processed in decreasing precedence order.

    Example:
        Run this function from the command line with the desired arguments, for example:

        ```bash
        sp_westpa_config --name my_config.cfg --save_dir /path/to/save --yaml config1.yaml config2.yaml
        ```

        This will generate a configuration file named "my_config.cfg" in the specified
        directory, updated with parameters from "config1.yaml" and "config2.yaml".
    """
    parser = argparse.ArgumentParser(description="Generate WESTPA configuration")
    parser.add_argument(
        "--name",
        type=str,
        help="File name including the file extension.",
        default="west.cfg",
    )
    parser.add_argument(
        "--save_dir", type=str, help="Directory to save rendered file.", default="."
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="+",
        help="Paths to YAML files to use for updating parameters in decreasing precedence.",
    )
    args = parser.parse_args()

    run_westpa_config(file_name=args.name, save_dir=args.save_dir, yaml_paths=args.yaml)
