"""Handle io of progress coordinate data."""

import os
from collections.abc import MutableSequence

from loguru import logger

from ..configs import SubpexConfig
from ..fop.io import write_fop


def _load_data(path_data: str) -> MutableSequence[float]:
    """Load and process the final progress coordinate value from a written file."""
    if os.path.exists(path_data):
        logger.info(f"Loading data at {path_data}")
        with open(path_data, "r", encoding="utf-8") as f:
            data_string = f.readlines()[-1].strip()

        data_split = data_string.split()
        if len(data_split) > 1:
            return [float(value) for value in data_split]
        return [float(data_split[0])]
    logger.debug(f"{path_data} does not exist")
    return []


def load_data(subpex_config: SubpexConfig, data_dir: str | None = None) -> SubpexConfig:
    """Initialize and load the last value of activated progress coordinate data.

    Args:
        subpex_config: Subpex context manager.
        data_dir: Directory to search for progress coordinate data files. If `None`,
            we will initialize the data with an empty list.

    Returns:
        The last value of data we could find.
    """
    if data_dir is not None:
        for aux_data in subpex_config.data.aux.get_active():
            _data = _load_data(os.path.join(data_dir, aux_data.file_name))
            if len(_data) > 0:
                aux_data.values.append(_data)

    return subpex_config


def write_data(
    subpex_config: SubpexConfig,
    data_dir: str | None = None,
) -> None:
    if data_dir is None:
        data_dir = ""
    for aux_data in subpex_config.data.aux.get_active():
        f_path = os.path.join(data_dir, aux_data.file_name)
        f_ext = f_path.split(".")[-1]

        if f_ext in ("xyz", "pdb"):
            if isinstance(aux_data.values[0], float):
                logger.warning(f"File extension for {f_path} is {f_ext}")
                logger.warning("However, data is not a field of points.")
                logger.warning("Not writing this data.")
                continue
            write_fop(
                aux_data.values,  # type: ignore
                f_path,
                data_dir=data_dir,
            )
        elif f_ext == "txt":
            if isinstance(aux_data.values[0], float):
                with open(f_path, "w", encoding="utf-8") as f:
                    for _value in aux_data.values:
                        f.write(str(_value) + "\n")
            elif isinstance(aux_data.values[0], list):
                with open(f_path, "w", encoding="utf-8") as f:
                    for _data_frame in aux_data.values:
                        line = ""
                        for _data_frame_pcoord in _data_frame:
                            line += "{:.4f}    ".format(_data_frame_pcoord)
                        f.write(line + "\n")
            else:
                logger.warning(f"{type(aux_data.values[0])} has not been accounted for")
                logger.warning(f"Skipping {f_path}")
        else:
            logger.warning(f"{f_ext} file extension is not supported")
            logger.warning(f"Skipping {f_path}")
