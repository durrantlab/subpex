import argparse
import os
import sys

from loguru import logger

from ...configs import SubpexConfig, WestpaConfig


def run_westpa_setup(
    subpex_config: SubpexConfig,
    westpa_config: WestpaConfig,
    write_dir: str = "",
    overwrite: bool = False,
) -> None:
    """Setup WESTPA simulation.

    Args:
        write_dir: Path to directory to initialize simulations.
    """
    if os.path.exists(write_dir):
        logger.info(f"Directory at {write_dir} exists")
        if overwrite:
            logger.info("Overwrite is `True`")
        else:
            logger.error("Overwrite is `False`")
            logger.error("Aborting")
            sys.exit(1)

    logger.info(f"Creating directory at {write_dir}")
    os.makedirs(write_dir, exist_ok=True)

    logger.info("Creating all necessary directories")
    logger.debug("   - Creating bstates")
    os.makedirs(os.path.join(write_dir, "bstates"))


def cli_run_westpa_setup():
    # Get arguments and load json settings file.
    parser = argparse.ArgumentParser(
        description="Initialize WESTPA simulation for SubPEx purposes."
    )
    parser.add_argument(
        "--write_dir",
        type=str,
        default="",
        help="Directory to write input files.",
    )
    parser.add_argument(
        "--yaml_subpex",
        type=str,
        nargs="+",
        help="Paths to YAML files to use for subpex config in decreasing precedence.",
    )
    parser.add_argument(
        "--yaml_westpa",
        type=str,
        nargs="+",
        help="Paths to YAML files to use for westpa config in decreasing precedence.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Remove and overwrite directory",
    )
    args = parser.parse_args()

    subpex_config = SubpexConfig()
    if args.yaml_subpex is None:
        args.yaml_subpex = []
    for yaml_path in reversed(args.yaml_subpex):
        subpex_config.from_yaml(yaml_path)

    westpa_config = WestpaConfig()
    if args.yaml_westpa is None:
        args.yaml_westpa = []
    for yaml_path in reversed(args.yaml_westpa):
        westpa_config.from_yaml(yaml_path)

    run_westpa_setup(subpex_config, westpa_config, args.write_dir)
