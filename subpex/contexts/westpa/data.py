from typing import Any

from collections.abc import Iterable, MutableMapping


class Data:
    """Data settings for WESTPA."""

    def __init__(self) -> None:
        self.west_data_file: str = "west.h5"
        """Name of the main HDF5 data storage file."""
        self.aux_compression_threshold: int = 1048576
        """Threshold in bytes for compressing auxiliary data."""
        self.iter_prec: int = 8
        """Length of the iteration index with zero-padding."""
        self.load_auxdata: bool = False
        """Load all datasets without an explicit load specification."""
        self.datasets: Iterable[MutableMapping[str, Any]] | None = None
        """List of all datasets to load."""
        self.data_refs: MutableMapping[str, str] = {
            "segment": r"$WEST_SIM_ROOT/traj_segs/{segment.n_iter:06d}/{segment.seg_id:06d}",
            "basis_state": r"$WEST_SIM_ROOT/bstates/{basis_state.auxref}",
            "initial_state": r"$WEST_SIM_ROOT/istates/{initial_state.iter_created}/{initial_state.state_id}.xml",
        }
        """References to various data paths."""
