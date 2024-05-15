from typing import Any

from collections.abc import Iterable


class ClusterContextManager:
    """Contexts for progress coordinate."""

    # pylint: disable-next=too-many-statements
    def __init__(
        self, yaml_paths: str | Iterable[str] | None = None, **kwargs: dict[str, Any]
    ) -> None:
        """
        Args:
            yaml_paths: Path(s) to YAML file(s) to load into the context.
        """
        self.clustering_engine: str = "cpptraj"
        """Clustering engine to use. `cpptraj` is the only option at the moment.
        """
        self.calculated_points: int = -1
        """Number of point to calculate per trajectory segment. If `-1`, it will
        calculate all.
        """
        self.n_clusters: int = 25
        """Number of independent clusters to identify."""
        self.min_n_clusters_gen_bin: int = 3
        """TODO:"""
        self.max_n_clusters_gen_bin: int = 25
        """TODO:"""
        self.clustering_region: str = "pocket"
        """Region to separate into clusters. This can be `pocket` or `backbone`."""
        self.clustering_method: str = "hierarchical"
        """TODO:"""
        self.cluster_generation: bool = True
        """TODO:"""
        self.cluster_bin: bool = False
        """TODO:"""
        self.plot_n_walkers: bool = False
        """TODO:"""
