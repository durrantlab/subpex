from typing import Any

import argparse
import os
import sys

from loguru import logger

from ...configs import SubpexConfig, WestpaConfig
from .scripts import write_env_script
from .utils import link_file


def link_initial_sims(
    subpex_config: SubpexConfig,
    westpa_config: WestpaConfig,
    write_dir: str = "",
    *args: Any,
    **kwargs: dict[str, Any],
) -> None:
    """WESTPA starts from a relaxed molecular simulation using the same MD engine.
    Instead of copying the preliminary simulation files, we create a soft link to
    the relevant files.
    """
    logger.info("Linking simulation files")
    link_file(
        path_original=subpex_config.sim.path_traj_relax,
        link_name="relax_traj",
        write_dir=os.path.join(write_dir, kwargs.get("common_files", "")),  # type: ignore
    )
    link_file(
        path_original=subpex_config.sim.path_topo,
        link_name="mol",
        write_dir=os.path.join(write_dir, kwargs.get("common_files", "")),  # type: ignore
    )


def create_westpa_dirs(write_dir: str = "", overwrite: bool = False) -> dict[str, str]:
    """
    Create necessary WESTPA directories.

    Args:
        write_dir (str): The directory where the WESTPA directories will be created. Defaults to an empty string.
        overwrite (bool): If True, existing directories will be overwritten. Defaults to False.

    Returns:
        dict[str, str]: A dictionary containing the paths of the created directories.
    """
    paths = {}
    logger.info("Creating all necessary WESTPA directories")
    logger.debug("   - Creating common_files")
    paths["common_files"] = os.path.join(write_dir, "common_files")
    os.makedirs(paths["common_files"], exist_ok=overwrite)
    logger.debug("   - Creating westpa_scripts")
    paths["westpa_scripts"] = os.path.join(write_dir, "westpa_scripts")
    os.makedirs(paths["westpa_scripts"], exist_ok=overwrite)
    logger.debug("   - Creating bstates")
    paths["bstates"] = os.path.join(write_dir, "bstates")
    os.makedirs(paths["bstates"], exist_ok=overwrite)
    return paths


def write_westpa_bash_scripts(
    subpex_config: SubpexConfig,
    westpa_config: WestpaConfig,
    write_dir: str = "",
    *args: Any,
    **kwargs: dict[str, Any],
) -> None:
    """
    Write WESTPA bash scripts based on the provided configurations.

    Args:
        subpex_config (SubpexConfig): The configuration for Subpex.
        westpa_config (WestpaConfig): The configuration for WESTPA.
        write_dir (str, optional): The directory where the scripts will be written. Defaults to "".
        *args: Additional positional arguments.
        **kwargs: Additional keyword arguments.

    Returns:
        None
    """
    write_env_script(westpa_config=westpa_config, write_dir=write_dir)


def run_westpa_setup(
    subpex_config: SubpexConfig,
    westpa_config: WestpaConfig,
    write_dir: str = "",
    overwrite: bool = False,
) -> None:
    """
    Run the setup process for WESTPA.

    Args:
        subpex_config (SubpexConfig): The configuration for the Subpex module.
        westpa_config (WestpaConfig): The configuration for the WESTPA module.
        write_dir (str, optional): The directory to write the setup files to. Defaults to "".
        overwrite (bool, optional): Whether to overwrite the existing directory if it exists. Defaults to False.
    """

    if os.path.exists(write_dir):
        logger.info(f"Directory at {write_dir} exists")
        if overwrite:
            logger.info("Overwrite is `True`")
        else:
            logger.error("Overwrite is `False`")
            logger.error("Aborting")
            sys.exit(1)

    # Create all directories
    logger.info(f"Creating directory at {write_dir} if absent")
    os.makedirs(write_dir, exist_ok=overwrite)

    path_dirs_westpa = create_westpa_dirs(write_dir=write_dir, overwrite=overwrite)

    # Link all necessary files
    link_initial_sims(
        subpex_config=subpex_config,
        westpa_config=westpa_config,
        write_dir=write_dir,
        **path_dirs_westpa,  # type: ignore
    )
    write_westpa_bash_scripts(
        subpex_config=subpex_config,
        westpa_config=westpa_config,
        write_dir=write_dir,
    )


def cli_run_westpa_setup():
    """
    Run the WESTPA setup for SubPEx purposes.

    This function initializes the WESTPA simulation for SubPEx purposes by parsing command line arguments,
    loading JSON settings file, and running the WESTPA setup with the provided configurations.

    Args:
        --write_dir (str): Directory to write input files.
        --yaml_subpex (List[str]): Paths to YAML files to use for subpex config in decreasing precedence.
        --yaml_westpa (List[str]): Paths to YAML files to use for westpa config in decreasing precedence.
        --overwrite (bool): Flag to overwrite files in write_dir.

    Returns:
        None
    """
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
        help="Overwrite files in write_dir",
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
