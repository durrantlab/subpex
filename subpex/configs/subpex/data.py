from collections.abc import MutableSequence

from pydantic import BaseModel, Field, computed_field


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

    @computed_field  # type: ignore[misc]
    @property
    def file_name(self) -> str:
        return self.label + "." + self.f_ext


class DataConfig(BaseModel):

    aux: MutableSequence[AuxData] = Field(
        default_factory=lambda: [
            AuxData(label="progress_coord"),
            AuxData(label="pocket_volume"),
            AuxData(label="pocket_rog"),
            AuxData(label="pocket_rmsd"),
            AuxData(label="pocket_fop"),
            AuxData(label="pocket_jd"),
            AuxData(label="backbone_rmsd"),
            AuxData(label="composite"),
        ]
    )

    write_dir: str = Field(default="")
    """Where to write any auxillary data that has `write` as `True`"""

    @property
    def progress_coord(self) -> MutableSequence[AuxData]:
        return [aux for aux in self.aux if aux.progress]

    """Auxillary data used in the progress coordinate"""

    def dict(self, **kwargs):
        d = super().dict(**kwargs)
        d["progress_coord"] = [aux.label for aux in self.progress_coord]
        return d

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
