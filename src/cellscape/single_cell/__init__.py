"""Single-cell visualization tools."""

from cellscape.single_cell.plots import (
    local_correlation_plot,
    umap_expr_with_category,
    umap_highlight,
)

__all__ = ["local_correlation_plot", "umap_expr_with_category", "umap_highlight"]
