from typing import Any, Generator

from collections.abc import Callable, MutableSequence

import MDAnalysis as mda
from pydantic import BaseModel, Field, computed_field

from ...pocket.compute import (
    get_pocket_jaccard_convenience,
    get_pocket_rmsd_convenience,
    get_pocket_rog_convenience,
    get_pocket_volume_convenience,
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
            result = self.compute_function(  # type: ignore
                atoms_frame, subpex_config, atoms_ref, *args, **kwargs
            )
            self.values.append(result)

    @computed_field  # type: ignore[misc]
    @property
    def file_name(self) -> str:
        return self.label + "." + self.f_ext


class ContainerAuxData(BaseModel):

    pocket_volume: AuxData = Field(
        default_factory=lambda: AuxData(
            label="pocket_volume", compute_function=get_pocket_volume_convenience  # type: ignore
        )
    )
    pocket_rog: AuxData = Field(
        default_factory=lambda: AuxData(
            label="pocket_rog", compute_function=get_pocket_rog_convenience  # type: ignore
        )
    )
    pocket_rmsd: AuxData = Field(
        default_factory=lambda: AuxData(
            label="pocket_rmsd", compute_function=get_pocket_rmsd_convenience  # type: ignore
        )
    )
    pocket_jd: AuxData = Field(
        default_factory=lambda: AuxData(
            label="pocket_jd", compute_function=get_pocket_jaccard_convenience  # type: ignore
        )
    )
    backbone_rmsd: AuxData = Field(
        default_factory=lambda: AuxData(
            label="backbone_rmsd", compute_function=get_backbone_rmsd  # type: ignore
        )
    )

    def get_all(self) -> Generator[AuxData, None, None]:
        """Yields each attribute of type AuxData."""
        for attribute in self.__annotations__:
            value = getattr(self, attribute)
            if isinstance(value, AuxData):
                yield value

    def get_active(self) -> Generator[AuxData, None, None]:
        """Yields each attribute of type AuxData."""
        for value in self.get_all():
            if value.active:
                yield value


class DataConfig(BaseModel):

    aux: ContainerAuxData = Field(default_factory=ContainerAuxData)

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
