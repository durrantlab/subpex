from pydantic import BaseModel, Field


class ClusteringConfig(BaseModel):

    clustering_engine: str = Field(default="cpptraj")
    """Clustering engine to use. `cpptraj` is the only option at the moment.
    """
    n_clusters: int = Field(default=25)
    """Number of independent clusters to identify."""
    min_n_clusters_gen_bin: int = Field(default=3)
    """TODO:"""
    max_n_clusters_gen_bin: int = Field(default=25)
    """TODO:"""
    clustering_region: str = Field(default="pocket")
    """Region to separate into clusters. This can be `pocket` or `backbone`."""
    clustering_method: str = Field(default="hierarchical")
    """TODO:"""
    cluster_generation: bool = Field(default=True)
    """TODO:"""
    cluster_bin: bool = Field(default=False)
    """TODO:"""
    plot_n_walkers: bool = Field(default=False)
    """TODO:"""
