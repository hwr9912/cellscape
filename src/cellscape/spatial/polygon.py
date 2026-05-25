"""Polygon and cell-boundary plotting for spatial omics data."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

from cellscape.core.plotting import apply_publication_defaults, finish_figure
from cellscape.core.validation import get_spatial_coordinates, require_obsm_key
from cellscape.spatial.plots import _get_values, _try_int


def cell_boundary_plot(
    adata: Any,
    color: str | None = None,
    *,
    border_key: str = "cell_border",
    spatial_key: str = "spatial",
    x_key: str = "x",
    y_key: str = "y",
    relative: bool = True,
    ax: plt.Axes | None = None,
    facecolor: str = "#DEDEDE",
    edgecolor: str = "none",
    linewidth: float = 0.1,
    alpha: float = 1.0,
    cmap: str = "viridis",
    palette: dict[Any, str] | None = None,
    show: bool = True,
    save: str | Path | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Plot cell polygons stored in ``adata.obsm[border_key]``.

    ``adata.obsm[border_key]`` must have shape ``(n_cells, n_vertices, 2)``.
    When ``relative=True``, each border is shifted by cell center coordinates
    from ``obs[x_key]``/``obs[y_key]`` or ``obsm[spatial_key]``.
    """
    apply_publication_defaults()
    require_obsm_key(adata, border_key)
    borders = np.asarray(adata.obsm[border_key])
    if borders.ndim != 3 or borders.shape[2] != 2:
        raise ValueError(
            f"adata.obsm[{border_key!r}] must have shape (n_cells, n_vertices, 2); "
            f"got {borders.shape}."
        )

    if relative:
        try:
            centers = get_spatial_coordinates(adata, x_key=x_key, y_key=y_key)
        except KeyError:
            centers = get_spatial_coordinates(adata, spatial_key=spatial_key)
        borders = borders + centers[:, None, :]

    patches = [Polygon(points, closed=True) for points in borders]
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.figure

    collection = PatchCollection(
        patches,
        linewidth=linewidth,
        edgecolor=edgecolor,
        alpha=alpha,
        match_original=False,
    )

    if color is None:
        collection.set_facecolor(facecolor)
    else:
        values = _get_values(adata, color)
        if pd.api.types.is_numeric_dtype(values):
            collection.set_array(values.to_numpy())
            collection.set_cmap(cmap)
            fig.colorbar(collection, ax=ax, fraction=0.046, pad=0.04)
        else:
            cats = pd.Series(values, dtype="category").cat.categories
            color_map = palette or {
                cat: plt.get_cmap("tab20")(i % 20) for i, cat in enumerate(cats)
            }
            collection.set_facecolor(values.map(color_map).fillna(facecolor))
            handles = [
                plt.Line2D([0], [0], marker="s", color=color_map[cat], linestyle="", label=str(cat))
                for cat in sorted(cats, key=_try_int)
            ]
            ax.legend(handles=handles, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")

    ax.add_collection(collection)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_axis_off()
    finish_figure(fig, save=save, show=show)
    return fig, ax
