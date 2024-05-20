"""Cluster field of points to identify pockets."""

from typing import Any

from collections.abc import MutableMapping, MutableSequence, Sequence

from sklearn.cluster import DBSCAN


def cluster_fop(
    fop: Sequence[Sequence[float]],
    clustering_model: type | None = None,
    clustering_kwargs: MutableMapping[str, Any] | None = None,
) -> Sequence[Sequence[float]]:
    """Takes the coordinates of the field of points (fop) and uses a clustering
    algorithm to trim the fop to the points. Keeps the cluster that contains
    the center of the pocket.

    Args:
        fop: coordinates of field of points.
        clustering_model: Scikit-learn clustering model class to use.
        clustering_kwargs: Keyword arguments for `clustering_model` initialization.

    Returns:
        field of points after clustering
    """
    # sometimes no points are present in the FOP.
    if len(fop) == 0:
        return []

    if clustering_model is None:
        clustering_model = DBSCAN
    if clustering_kwargs is None:
        if isinstance(clustering_model, DBSCAN):
            clustering_kwargs = {"eps": 0.5, "min_samples": 2}
        else:
            clustering_kwargs = {}
    db = clustering_model(**clustering_kwargs).fit(fop)
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    if n_clusters_ == 0:
        # In some cases, even when fop has a few points, there are no clusters.
        # In this case, just return an empty list.
        return []

    longest: str = ""
    clusters: dict[str, MutableSequence[Sequence[float]]] = {}
    for i in range(n_clusters_):
        clusters["Cluster_" + str(i + 1)] = []

    for i, coord in enumerate(fop):
        if labels[i] != -1:
            clusters["Cluster_" + str(labels[i] + 1)].append(coord)

    for cluster_key, cluster in clusters.items():
        if longest == "" or len(cluster) > len(clusters[longest]):
            longest = cluster_key

    return clusters[longest]
