from typing import Any

from collections.abc import Iterable, MutableSequence


class DataContextManager:
    """Contexts for progress coordinate."""

    def __init__(
        self, yaml_paths: str | Iterable[str] | None = None, **kwargs: dict[str, Any]
    ) -> None:
        """
        Args:
            yaml_paths: Path(s) to YAML file(s) to load into the context.
        """
        self.progress_coord: MutableSequence[str] = ["composite"]
        """Progress coordinates to calculate.

        -   `jd` for Jaccard distance,
        -   `prmsd` for pocket heavy atoms RMSD,
        -   `bb` for backbone RMSD,
        -   `composite` for composite RMSD.
        """
        self.clustering_engine: str = "cpptraj"
        """Clustering engine to use. `cpptraj` is the only option at the moment.
        """
        self.calculated_points: int = -1
        """Number of point to calculate per trajectory segment. If `-1`, it will
        calculate all.
        """
        self.ref_fop_write: str | None = "fop-ref.xyz"
        """Path to write the reference field of points."""
        self.ref_pdb_path: str | None = None
        """Path to a PDB file to use as a reference to compute the progress coordinate.
        """

        self.aux_progress_coord: bool = True
        """Flag for computing the progress coordinate."""
        self.aux_pocket_vol: bool = True
        """Flag for computing the pocket volume data."""
        self.aux_pocket_rog: bool = True
        """Flag for computing the pocket radius of gyration."""
        self.aux_pocket_rmsd: bool = True
        """Flag for computing the pocket RMSD."""
        self.aux_pocket_fop: bool = True
        """Flag for computing the the pocket field of points data."""
        self.aux_pocket_jd: bool = True
        """Flag for computing the the Jaccard distance of the pocket."""
        self.aux_bb_rmsd: bool = True
        """Flag for computing the backbone RMSD."""
        self.aux_composite: bool = True
        """Flag for computing the pocket radius of gyration."""

        self.key_progress_coord: str = "progress_coords"
        """Key for the progress coordinate."""
        self.key_pocket_vol: str = "pocket_vol"
        """Key for pocket volume data."""
        self.key_pocket_rog: str = "pocket_rog"
        """Key for pocket radius of gyration."""
        self.key_pocket_rmsd: str = "pocket_rmsd"
        """Key for pocket RMSD."""
        self.key_pocket_fop: str = "pocket_fop"
        """Key for the pocket field of points data."""
        self.key_pocket_jd: str = "pocket_jd"
        """Key for the Jaccard distance of the pocket."""
        self.key_bb_rmsd: str = "bb_rmsd"
        """Key for backbone RMSD."""
        self.key_composite: str = "composite"
        """Key for pocket radius of gyration."""

        self.fname_progress_coord: str = "progress_coord.txt"
        """File name of the progress coordinate."""
        self.fname_pocket_vol: str = "pocket_vol.txt"
        """File name of pocket volume data."""
        self.fname_pocket_rog: str = "pocket_rog.txt"
        """File name of pocket radius of gyration."""
        self.fname_pocket_rmsd: str = "pocket_rmsd.txt"
        """File name of pocket RMSD."""
        self.fname_pocket_fop: str = "pocket_fop.txt"
        """File name of the pocket field of points data."""
        self.fname_pocket_jd: str = "pocket_jd.txt"
        """File name of the Jaccard distance of the pocket."""
        self.fname_bb_rmsd: str = "bb_rmsd.txt"
        """File name of backbone RMSD."""
        self.fname_composite: str = "composite.txt"
        """File name of pocket radius of gyration."""

    @property
    def data_keys(self) -> MutableSequence[str]:
        """Keys for possible data."""
        keys = [attr for attr in vars(self).keys() if attr.startswith("key_")]
        return keys

    @property
    def data_file_names(self) -> MutableSequence[str]:
        """File names for possible data. These include the file extension."""
        fnames = [attr for attr in vars(self).keys() if attr.startswith("fname_")]
        return fnames

    @property
    def data_activated(self) -> MutableSequence[str]:
        """File names for possible data. These include the file extension."""
        activated = [attr for attr in vars(self).keys() if attr.startswith("aux_")]
        return activated
