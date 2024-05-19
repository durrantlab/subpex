"""Handle io of progress coordinate data."""

import os
from collections.abc import MutableMapping, MutableSequence

from loguru import logger

from ..configs import SubpexConfig
from ..fop.io import write_fop


def _load_data(data_dir: str, fname_data: str) -> float | MutableSequence[float]:
    """Load and process the final progress coordinate value from a written file."""
    path_data = os.path.join(data_dir, fname_data)
    if os.path.exists(path_data):
        logger.info(f"Loading data at {path_data}")
        with open(path_data, "r", encoding="utf-8") as f:
            data_string = f.readlines()[-1].strip()

        data_split = data_string.split()
        if len(data_split) > 1:
            return [float(value) for value in data_split]
        return float(data_split[0])
    logger.debug(f"{path_data} does not exist")
    return []


def initialize_data(
    subpex_config: SubpexConfig, data_dir: str | None = None
) -> MutableMapping[str, MutableSequence[float | MutableSequence[float]]]:
    """Initialize and load the last value of activated progress coordinate data.

    Args:
        subpex_config: Subpex context manager.
        data_dir: Directory to search for progress coordinate data files. If `None`,
            we will initialize the data with an empty list.

    Returns:
        The last value of data we could find.
    """
    data: MutableMapping[str, MutableSequence[float | MutableSequence[float]]] = {}

    for activated, key, fname in zip(
        sorted(subpex_config.data_activated),
        sorted(subpex_config.data_keys),
        sorted(subpex_config.data_file_names),
    ):
        if not activated:
            continue
        _data = []
        if data_dir is not None:
            _data.append(_load_data(data_dir, getattr(subpex_config, fname)))
        data[getattr(subpex_config, key)] = _data

    return data


def write_data(
    data: MutableMapping[str, MutableSequence[float | MutableSequence[float]]],
    subpex_config: SubpexConfig,
    data_dir: str | None = None,
) -> None:
    if data_dir is None:
        data_dir = ""
    for activated, key, fname in zip(
        sorted(subpex_config.data_activated),
        sorted(subpex_config.data_keys),
        sorted(subpex_config.data_file_names),
    ):
        if not activated:
            continue
        f_path = os.path.join(data_dir, fname)
        f_ext = fname.split(".")[-1]

        _data = data[getattr(subpex_config, key)]
        if f_ext in ("xyz", "pdb"):
            if isinstance(_data[0], float):
                logger.warning(f"File extension for {key} is {f_ext}")
                logger.warning("However, data is not a field of points.")
                logger.warning("Not writing this data.")
                continue
            write_fop(
                data[getattr(subpex_config, key)],
                f_path,
                subpex_config=subpex_config,
                data_dir=data_dir,
            )
        elif f_ext == "txt":
            if isinstance(_data[0], float):
                with open(f_path, "w", encoding="utf-8") as f:
                    for _value in _data:
                        f.write(str(_value) + "\n")
            elif isinstance(_data[0], list):
                with open(f_path, "w", encoding="utf-8") as f:
                    for _data_frame in _data:
                        line = ""
                        for _data_frame_pcoord in _data_frame:
                            line += "{:.4f}    ".format(_data_frame_pcoord)
                        f.write(line + "\n")
            else:
                logger.warning(f"{type(_data[0])} has not been accounted for")
                logger.warning(f"Skipping {key}")
        else:
            logger.warning(f"{f_ext} file extension is not supported")
            logger.warning(f"Skipping {key}")
