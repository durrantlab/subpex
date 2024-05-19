"""Assist with rendering templates"""

import argparse
import os
import re
import sys
from collections.abc import Sequence

from jinja2 import Template

from .. import DIR_TEMPLATE
from ..configs import SubpexConfig, WestpaConfig

SLURM_TEMPLATE_PATH = os.path.join(DIR_TEMPLATE, "workload_managers/job.sbatch")


def load_template(file_path: str) -> str:
    with open(file_path, "r", encoding="utf-8") as f:
        return f.read()


def clean_render(render_string: str) -> str:
    # Remove leading and trailing whitespace
    text = render_string.strip()
    # Replace multiple consecutive empty lines with a single empty line
    text = re.sub(r"\n\s*\n+", "\n\n", text)
    # Ensure there is exactly one trailing newline
    return text + "\n"


def slurm(
    subpex_config: SubpexConfig,
    file_name: str | None = None,
    save_dir: str | None = None,
) -> Sequence[str]:
    """Render slurm sbatch script.

    Args:
        subpex_config: SubPEx context manager.
        file_name: Name for the file including the extension.
        save_dir: Directory to save rendered file.
    """
    if file_name is None:
        file_name = "job.sbatch"

    template = Template(load_template(SLURM_TEMPLATE_PATH))

    rendered_file = template.render(**subpex_config.get())
    rendered_file = clean_render(rendered_file)

    if save_dir is not None:
        save_path = os.path.join(save_dir, file_name)
        with open(save_path, "w", encoding="utf-8") as f:
            f.write(rendered_file)

    return rendered_file


def west_cfg(
    westpa_cm: WestpaConfig,
    file_name: str | None = None,
    save_dir: str | None = None,
) -> None:
    """Writes `west.cfg` file for WESTPA simulations."""
    westpa_cm.to_yaml(save_dir=save_dir)


def cli_renderer():
    parser = argparse.ArgumentParser(description="Render templates provided by SubPEx.")
    parser.add_argument(
        "function_name",
        type=str,
        help="Rendering function to use in `subpex.setup.render`.",
    )
    parser.add_argument(
        "--name", type=str, help="File name including the file extension."
    )
    parser.add_argument("--save_dir", type=str, help="Directory to save rendered file.")
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="+",
        help="Paths to YAML files to use for SubpexContext in decreasing precedence.",
    )
    args = parser.parse_args()

    if args.function_name == "west_cfg":
        cm = WestpaConfig()
    elif args.function_name in ["slurm"]:
        cm = SubpexConfig()
    else:
        raise ValueError(f"Could not find function with name {args.function_name}")
    if args.yaml is None:
        args.yaml = []
    for yaml_path in reversed(args.yaml):
        cm.from_yaml(yaml_path)

    if args.save_dir is None:
        args.save_dir = ""

    current_module = sys.modules[__name__]
    renderer = getattr(current_module, args.function_name, None)
    if callable(renderer):
        return renderer(cm, args.name, args.save_dir)
    else:
        raise ValueError(f"No callable function named '{args.function_name}' found.")
