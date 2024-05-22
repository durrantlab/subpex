from typing import Any

from collections.abc import Iterable, MutableMapping

from pydantic import BaseModel, Field


class DatasetsConfig(BaseModel):
    name: str
    h5path: str | None = Field(default=None)
    store: bool = Field(default=True)
    load: bool = Field(default=False)
    dtype: str | None = Field(default=None)
    scaleoffset: int | None = Field(default=None)
    compression: int | None = Field(default=None)
    chunks: int | None = Field(default=None)


class DataRefs(BaseModel):
    segment: str | None = Field(
        default=r"$WEST_SIM_ROOT/traj_segs/{segment.n_iter:06d}/{segment.seg_id:06d}"
    )
    basis_state: str | None = Field(
        default=r"$WEST_SIM_ROOT/bstates/{basis_state.auxref}"
    )
    initial_state: str | None = Field(
        default=r"$WEST_SIM_ROOT/istates/{initial_state.iter_created}/{initial_state.state_id}.xml"
    )


class DataConfig(BaseModel):

    west_data_file: str = Field(default="west.h5")
    """The name of the main HDF5 data storage file for the WESTPA simulation."""
    aux_compression_threshold: int = Field(default=1048576)
    """The threshold in bytes for compressing the auxiliary data in a dataset on an
    iteration-by-iteration basis.
    """
    iter_prec: int = Field(default=8)
    """The length of the iteration index with zero-padding. For the default value,
    iteration 1 would be specified as iter_00000001.
    """
    load_auxdata: bool = Field(default=False)
    """Load all datasets without an explicit load specification."""
    datasets: Iterable[MutableMapping[str, Any]] | None = Field(default=None)
    """List of all datasets to load."""
    data_refs: DataRefs = Field(default_factory=DataRefs)
    """References to various data paths."""
