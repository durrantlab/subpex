from collections.abc import Sequence

from pydantic import BaseModel, Field


class PocketConfig(BaseModel):
    """Contexts for progress coordinate."""

    selection_mda_str: str | None = Field(default=None)
    """MDAnalysis atom selection string for pocket.
    """
    selection_dist: float | None = Field(default=None)
    """Define the pocket as all atoms within a radius of the center.
    """
    selection_append: str | None = "and protein and (not name H*)"
    """String to append to pocket selection if desired."""
    resids: Sequence[int] | None = Field(default=None)
    """Residue IDs to form the pocket.
    """
    selection_str: str | None = Field(default=None)
    """Final pocket selection string."""
    pocket_center: Sequence[float] | None = Field(default=None)
    """The pocket center XYZ coordinates."""
    radius: float = Field(default=10.0)
    """Pocket radius to consider."""
    resolution: float = Field(default=0.5)
    """Resolution for the FOP points."""
    selection_write: str | None = Field(default=None)
    """Path to write a pocket selection string text file. If `None`, then no file
    is written.
    """
    ref_fop_write: str | None = "fop_ref.xyz"
    """Path to write the reference field of points."""
    ref_pdb_path: str | None = None
    """Path to a PDB file to use as a reference to compute the progress coordinate.
    """
