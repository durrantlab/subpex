"""Handle io of progress coordinate data."""

import os
from collections.abc import MutableMapping, MutableSequence

from loguru import logger

from ..contexts.subpex import SubpexContextManager


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


def initialize_pcoord_data(
    subpex_cm: SubpexContextManager, data_dir: str | None = None
) -> MutableMapping[str, MutableSequence[float] | float]:
    """Initialize and load the last value of activated progress coordinate data.

    Args:
        subpex_cm: Subpex context manager.
        data_dir: Directory to search for progress coordinate data files. If `None`,
            we will initialize the data with an empty list.

    Returns:
        The last value of data we could find.
    """
    data = {}

    for activated, key, fname in zip(
        subpex_cm.data_activated, subpex_cm.data_keys, subpex_cm.data_file_names
    ):
        if not activated:
            continue
        if data_dir is not None:
            _data = _load_data(data_dir, getattr(subpex_cm, fname))
        else:
            _data = []
        data[getattr(subpex_cm, key)] = _data

    return data
