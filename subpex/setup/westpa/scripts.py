"""Manages setting up scripts for WESTPA simulations."""

import os

from loguru import logger

from .configs import WestpaConfig


def write_env_script(westpa_config: WestpaConfig, write_dir: str = "") -> None:
    """
    Write the environment script (`env.sh`) based on the provided WestpaConfig object.

    Args:
        westpa_config: The WestpaConfig object containing the environment variables.
        write_dir: The directory where the env.sh file will be written.
    """
    logger.debug("Writing env.sh")
    lines_env = []
    for key, value in westpa_config.env.dict().items():
        if value is None:
            logger.trace(f"Skipping {key} from env.sh because it is `None`.")
            continue
        logger.trace(f"Writing {key}={value} to env.sh")
        lines_env.append(f"export {key}={value}")
    with open(os.path.join(write_dir, "env.sh"), "w") as f:
        f.write("\n".join(lines_env))
    os.chmod(os.path.join(write_dir, "env.sh"), 0o777)
