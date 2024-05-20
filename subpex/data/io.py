"""Handle IO of auxillary data."""

import os
from collections.abc import MutableSequence

from loguru import logger

from ..configs import SubpexConfig
from ..fop.io import write_fop


def _load_data(path_data: str) -> MutableSequence[float]:
    """Attempts to load the last auxillary data from a data file, if it exists.

    This function reads a specified file and retrieves the last value, which is assumed
    to be stored on the last line of the file. If the path does not exist, it will
    return an empty list.

    Args:
        path_data: The path to the file containing progress coordinate data.
    """
    if os.path.exists(path_data):
        logger.debug(f"Loading data at {path_data}")
        with open(path_data, "r", encoding="utf-8") as f:
            data_string = f.readlines()[-1].strip()

        data_split = data_string.split()
        if len(data_split) > 1:
            return [float(value) for value in data_split]
        return [float(data_split[0])]
    logger.debug(f"{path_data} does not exist")
    return []


def load_data(subpex_config: SubpexConfig, data_dir: str | None = None) -> SubpexConfig:
    """Populate `subpex_config` auxillary data with the last value found in `data_dir`.

    Calls

    Args:
        subpex_config: The Subpex configuration.
        data_dir: Directory to search for auxillary data files. If `None`, then no
            data is appended to `subpex_config`.

    Returns:
        The updated SubpexConfig object with the last value of data found appended to
        the respective auxiliary data.
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
    """Writes activated auxiliary data from the `subpex_config` to files in the
    specified directory.

    This function iterates through the active auxiliary data in the provided
    `subpex_config` and writes them to files in the `data_dir`. The format of the output
    files is determined by their extensions:

    - `.xyz` and `.pdb`: The data is expected to be a field of points and is written
      using the `write_fop` function.
    - `.txt`: The data is written as plain text. If the data consists of floats, each
      value is written on a new line. If the data is iterable (e.g., a list of lists),
      each sublist is written on a new line, with values formatted to four decimal
      places and separated by spaces.

    Args:
        subpex_config: The Subpex configuration object containing auxiliary data to be
            written to files.
        data_dir: The directory where the data files will be written. If `None`, an
            empty string is used as the directory.

    Warning:
        - If the file extension is not supported (`.xyz`, `.pdb`, or `.txt`), the
          function logs a warning and skips writing the data.
        - If the data type does not match the expected format for a given file
          extension, the function logs a warning and skips writing the data.
        - If the data type is not accounted for, the function logs a warning and skips
          writing the data.

    Example:
        ```python
        subpex_config = SubpexConfig(...)
        write_data(subpex_config, data_dir="/path/to/data")
        ```
    """
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
            elif hasattr(aux_data.values[0], "__iter__"):
                with open(f_path, "w", encoding="utf-8") as f:
                    for _data_frame in aux_data.values:
                        line = ""
                        for _data_frame_pcoord in _data_frame:  # type: ignore
                            line += "{:.4f}    ".format(_data_frame_pcoord)
                        f.write(line + "\n")
            else:
                logger.warning(f"{type(aux_data.values[0])} has not been accounted for")
                logger.warning(f"Skipping {f_path}")
        else:
            logger.warning(f"{f_ext} file extension is not supported")
            logger.warning(f"Skipping {f_path}")
