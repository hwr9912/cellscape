"""Visualization tools for single-cell and spatial omics analysis."""

from importlib.metadata import PackageNotFoundError, version

from cellscape.spatial.annotation import project_labelme_masks_to_obs

try:
    __version__ = version("cellscape")
except PackageNotFoundError:
    __version__ = "0.0.0"

__all__ = [
    "__version__",
    "LabelmeMaskResult",
    "batch_labelme_to_masks",
    "cell_boundary_plot",
    "cellscape_continuous_cmap",
    "highlight_and_expression_grid",
    "highlight_clusters_panels",
    "labelme_to_mask",
    "labelme_to_masks",
    "local_correlation_plot",
    "project_labelme_masks_to_obs",
    "spatial_expression_panels",
    "spatial_scatter",
    "umap_expr_with_category",
    "umap_highlight",
]


def __getattr__(name: str):
    """Lazily expose plotting helpers without drawing or loading data on import."""
    if name in {
        "cell_boundary_plot",
        "highlight_and_expression_grid",
        "highlight_clusters_panels",
        "spatial_expression_panels",
        "spatial_scatter",
    }:
        from cellscape import spatial

        return getattr(spatial, name)
    if name in {
        "LabelmeMaskResult",
        "batch_labelme_to_masks",
        "labelme_to_mask",
        "labelme_to_masks",
    }:
        from cellscape import datasets

        return getattr(datasets, name)
    if name in {
        "local_correlation_plot",
        "umap_expr_with_category",
        "umap_highlight",
    }:
        from cellscape import single_cell

        return getattr(single_cell, name)
    if name == "cellscape_continuous_cmap":
        from cellscape.styles import cellscape_continuous_cmap

        return cellscape_continuous_cmap
    raise AttributeError(f"module 'cellscape' has no attribute {name!r}")
