from typing import Any

from collections.abc import Callable, MutableSequence

import MDAnalysis as mda
from pydantic import BaseModel, Field, computed_field

from ...fop.compute import get_fop_volume_convenience
from ...pocket.compute import (
    get_pocket_jaccard_convenience,
    get_pocket_rmsd_convenience,
    get_pocket_rog_convenience,
)
from ...protein.compute import get_backbone_rmsd


class AuxData(BaseModel):
    label: str
    """Label used for any data structures or file names."""
    active: bool = Field(default=False)
    """If this data should be calculated and monitored during simulations or analyzed."""
    progress: bool = Field(default=False)
    """If this auxillary data should be used for a progress coordinate"""
    values: MutableSequence[float | MutableSequence[float]] = Field(
        default=[], exclude=True
    )
    """Data values."""
    f_ext: str = Field(default="txt")
    """File extension if saving."""
    write: bool = Field(default=False)
    """Flag to write file containing data."""
    compute_function: Callable[[Any], float | MutableSequence[float]] = Field(
        default=None, exclude=True
    )
    """Function to compute the descriptor value."""

    def compute_frame(  # type: ignore
        self,
        atoms_frame: mda.AtomGroup,
        subpex_config,
        atoms_ref: mda.AtomGroup | None = None,
        *args: Any,
        **kwargs: Any
    ) -> None:
        """Will use `compute_function` to compute and then append the value to
        `values`.
        """
        if self.compute_function is not None:
            result = self.compute_function(
                atoms_frame, subpex_config, atoms_ref, *args, **kwargs
            )
            self.values.append(result)

    @computed_field  # type: ignore[misc]
    @property
    def file_name(self) -> str:
        return self.label + "." + self.f_ext


class DataConfig(BaseModel):

    aux: MutableSequence[AuxData] = Field(
        default_factory=lambda: [
            AuxData(label="pocket_volume", compute_function=get_fop_volume_convenience),
            AuxData(label="pocket_rog", compute_function=get_pocket_rog_convenience),
            AuxData(label="pocket_rmsd", compute_function=get_pocket_rmsd_convenience),
            AuxData(label="pocket_jd", compute_function=get_pocket_jaccard_convenience),
            AuxData(label="backbone_rmsd", compute_function=get_backbone_rmsd),
        ]
    )

    write_dir: str = Field(default="")
    """Where to write any auxillary data that has `write` as `True`"""

    composite_sigma: float | None = Field(default=None)
    """Sigma to use for composite progress coordinate. If `None`, this defaults
    to

    ```python
    sigma = (
        1
        - len(reference.select_atoms(selection_alignment))
        / len(reference.select_atoms("backbone"))
    ) / 2
    ```
    """
