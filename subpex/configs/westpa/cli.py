import argparse
import os
from collections.abc import Sequence

from .main import WestpaConfig


def run_westpa_config(
    file_name: str = "west.cfg",
    save_dir: str | None = None,
    yaml_paths: Sequence[str] | None = None,
) -> None:
    # Initialize default System instance
    westpa_config = WestpaConfig()

    # Update parameters from YAML files if provided
    if yaml_paths:
        westpa_config.from_yaml(yaml_paths)

    # Create the output file path
    if save_dir is None:
        save_dir = ""
    output_path = os.path.join(save_dir, file_name)

    # Write the YAML file
    westpa_config.to_yaml(output_path)


def cli_westpa_config():
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

    run_westpa_config(save_dir=args.save_dir, yaml_paths=args.yaml)
