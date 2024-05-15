from typing import Any

from collections.abc import Iterable, Sequence


class PocketContextManager:
    """Contexts for progress coordinate."""

    # pylint: disable-next=too-many-statements
    def __init__(
        self, yaml_paths: str | Iterable[str] | None = None, **kwargs: dict[str, Any]
    ) -> None:
        """
        Args:
            yaml_paths: Path(s) to YAML file(s) to load into the context.
        """
        self.pocket_selection_mda: str | None = None
        """MDAnalysis atom selection string for pocket.
        """
        self.pocket_selection_resids: Sequence[int] | None = None
        """Residue IDs to form the pocket.
        """
        self.pocket_selection_dist: float | None = None
        """Define the pocket as all atoms within a radius of the center.
        """
        self.pocket_selection_append: str | None = "and protein and (not name H*)"
        """String to append to pocket selection if desired."""
        self.pocket_selection_str: str | None = None
        """Final pocket selection string."""
        self.pocket_center: Sequence[float] | None = None
        """The pocket center XYZ coordinates."""
        self.pocket_radius: float = 10.0
        """Pocket radius to consider."""
        self.pocket_resolution: float = 0.5
        """Resolution for the FOP points."""
        self.pocket_selection_write: str | None = None
        """Path to write a pocket selection string text file. If `None`, then no file
        is written.
        """
        self.ref_fop_write: str | None = "fop_ref.xyz"
        """Path to write the reference field of points."""
        self.ref_pdb_path: str | None = None
        """Path to a PDB file to use as a reference to compute the progress coordinate.
        """
