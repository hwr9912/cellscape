"""Spatial transcriptomics visualization tools."""

from cellscape.spatial.annotation import project_labelme_masks_to_obs

__all__ = [
    "cell_boundary_plot",
    "highlight_and_expression_grid",
    "highlight_clusters_panels",
    "project_labelme_masks_to_obs",
    "spatial_expression_panels",
    "spatial_scatter",
]


def __getattr__(name: str):
    """按需导入 spatial 子模块, 避免导入非绘图函数时加载完整绘图依赖"""
    if name == "cell_boundary_plot":
        from cellscape.spatial.polygon import cell_boundary_plot

        return cell_boundary_plot
    if name in {
        "highlight_and_expression_grid",
        "highlight_clusters_panels",
        "spatial_expression_panels",
        "spatial_scatter",
    }:
        from cellscape.spatial import plots

        return getattr(plots, name)
    raise AttributeError(f"module 'cellscape.spatial' has no attribute {name!r}")
