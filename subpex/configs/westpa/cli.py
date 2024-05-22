import argparse
import os
from collections.abc import Sequence

from loguru import logger

from .main import WestpaConfig


def run_westpa_config(
    westpa_config: WestpaConfig | None = None,
    file_name: str = "west.cfg",
    save_dir: str | None = None,
    yaml_paths: Sequence[str] | None = None,
) -> None:
    """
    Run the WestpaConfig configuration process.

    Args:
        westpa_config (WestpaConfig | None): An instance of WestpaConfig. If `None`,
            a default instance will be initialized.
        file_name (str): The name of the output file. Default is 'west.cfg'.
        save_dir (str | None): The directory where the output file will be saved.
            If `None`, the current directory will be used.
        yaml_paths (Sequence[str] | None): A sequence of YAML file paths to update the
            configuration parameters from.

    Returns:
        None: This function does not return anything.

    """
    if westpa_config is None:
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
