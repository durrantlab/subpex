from typing import Any

import os
from collections.abc import Iterable, MutableMapping

import yaml


class System:
    """System settings for WESTPA."""

    def __init__(self) -> None:
        self.driver: str = "adaptive.System"
        """The driver parameter must be set to a subclass of `WESTSystem`."""
        self.module_path: str = "$WEST_SIM_ROOT/adaptive_binning/"
        """Path where the class is defined."""


class Propagation:
    """Propagation settings for WESTPA."""

    def __init__(self) -> None:
        self.gen_istates: bool = True
        """Generate initial states from basis states."""
        self.block_size: int = 1
        """Number of segments passed to a worker at a time."""
        self.save_transition_matrices: bool = True
        """Save transition matrices."""
        self.max_run_wallclock: str | None = None
        """Max wallclock time for a WESTPA run."""
        self.max_total_iterations: int | None = None
        """Max number of iterations to run."""


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


class Executable:
    """Executable settings for WESTPA."""

    def __init__(self):
        self.pre_iteration: MutableMapping[str, Any] = {
            "enabled": False,
            "executable": "$WEST_SIM_ROOT/westpa_scripts/pre_iter.sh",
            "stderr": "$WEST_SIM_ROOT/job_logs/pre_iter.log",
        }
        """Pre-iteration executable settings."""
        self.post_iteration: MutableMapping[str, Any] = {
            "enabled": True,
            "executable": "$WEST_SIM_ROOT/westpa_scripts/post_iter.sh",
            "stderr": "$WEST_SIM_ROOT/job_logs/post_iter.log",
        }
        """Post-iteration executable settings."""


def to_dict_helper(obj: Any) -> Any:
    if hasattr(obj, "__dict__"):
        return {k: to_dict_helper(v) for k, v in obj.__dict__.items()}
    elif isinstance(obj, list):
        return [to_dict_helper(v) for v in obj]
    elif isinstance(obj, dict):
        return {k: to_dict_helper(v) for k, v in obj.items()}
    else:
        return obj


class WestpaContextManager:

    def __init__(self, yaml_paths: str | Iterable[str] | None = None, **kwargs: Any):
        self.system = System()
        self.propagation = Propagation()
        self.data = Data()
        self.executable = Executable()

        if yaml_paths:
            self.from_yaml(yaml_paths)
        self.update_config(**kwargs)

    def update_config(self, **kwargs: MutableMapping[str, Any]) -> None:
        for key, value in kwargs.items():
            keys = key.split(".")
            current_level = self
            for k in keys[:-1]:
                current_level = getattr(current_level, k)
            setattr(current_level, keys[-1], value)

    def from_yaml(self, paths: str | Iterable[str]) -> None:
        if isinstance(paths, str):
            paths = [paths]
        for path in paths:
            with open(path, "r") as f:
                config = yaml.safe_load(f)
                self.update_config(**config)

    def get(self):
        context = to_dict_helper(self)
        return context

    def to_yaml(self, save_dir: str | None = None) -> None:
        if save_dir is None:
            save_dir = ""
        yaml_path = os.path.join(save_dir, "west.cfg")
        with open(yaml_path, "w", encoding="utf-8") as f:
            yaml.dump({"west": self.get()}, stream=f, default_flow_style=False)
