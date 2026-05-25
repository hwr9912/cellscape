"""High-level plotting functions for spatial omics analysis."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cellscape.core.plotting import (
    apply_publication_defaults,
    as_axes_array,
    finish_figure,
)
from cellscape.core.validation import get_spatial_coordinates, require_color_key


def _as_list(value: str | Sequence[str]) -> list[str]:
    if isinstance(value, str):
        return [value]
    return list(value)


def _try_int(value: Any) -> Any:
    try:
        return int(value)
    except (TypeError, ValueError):
        return value


def _get_values(adata: Any, color: str, *, layer: str | None = None) -> pd.Series:
    if color in adata.obs:
        return adata.obs[color]
    if color not in adata.var_names:
        raise KeyError(f"{color!r} was not found in adata.obs or adata.var_names.")

    gene_idx = list(adata.var_names).index(color)
    matrix = adata.layers[layer] if layer is not None else adata.X
    values = matrix[:, gene_idx]
    if hasattr(values, "toarray"):
        values = values.toarray()
    return pd.Series(np.asarray(values).reshape(-1), index=adata.obs_names, name=color)


def spatial_scatter(
    adata: Any,
    color: str,
    *,
    spatial_key: str = "spatial",
    x_key: str | None = None,
    y_key: str | None = None,
    layer: str | None = None,
    ax: plt.Axes | None = None,
    size: float = 8,
    alpha: float = 1.0,
    cmap: str = "viridis",
    palette: dict[Any, str] | Sequence[str] | None = None,
    na_color: str = "#DEDEDE",
    title: str | None = None,
    frameon: bool = False,
    legend: bool = True,
    show: bool = True,
    save: str | Path | None = None,
    **scatter_kwargs: Any,
) -> tuple[plt.Figure, plt.Axes]:
    """Plot spatial coordinates colored by an obs column or gene expression."""
    apply_publication_defaults()
    require_color_key(adata, color)
    coords = get_spatial_coordinates(
        adata,
        spatial_key=spatial_key,
        x_key=x_key,
        y_key=y_key,
    )
    values = _get_values(adata, color, layer=layer)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    else:
        fig = ax.figure

    if pd.api.types.is_numeric_dtype(values):
        scatter = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=values.to_numpy(),
            s=size,
            alpha=alpha,
            cmap=cmap,
            **scatter_kwargs,
        )
        if legend:
            fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    else:
        cats = pd.Series(values, dtype="category").cat.categories
        if isinstance(palette, dict):
            color_map = {cat: palette.get(cat, na_color) for cat in cats}
        elif palette is not None:
            color_map = dict(zip(cats, palette))
        elif f"{color}_colors" in adata.uns and len(adata.uns[f"{color}_colors"]) == len(cats):
            color_map = dict(zip(cats, adata.uns[f"{color}_colors"]))
        else:
            cmap_obj = plt.get_cmap("tab20")
            color_map = {cat: cmap_obj(i % cmap_obj.N) for i, cat in enumerate(cats)}

        point_colors = values.map(color_map).fillna(na_color)
        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=point_colors,
            s=size,
            alpha=alpha,
            **scatter_kwargs,
        )
        if legend:
            handles = [
                mpatches.Patch(color=color_map[cat], label=str(cat))
                for cat in sorted(cats, key=_try_int)
            ]
            ax.legend(handles=handles, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")

    ax.set_title(title or color)
    ax.set_aspect("equal")
    if not frameon:
        ax.set_axis_off()
    finish_figure(fig, save=save, show=show)
    return fig, ax


def spatial_expression_panels(
    adata: Any,
    gene: str,
    panels: Sequence[str],
    *,
    library_key: str = "library",
    spatial_key: str = "spatial",
    layer: str | None = "lognorm",
    size: float = 50,
    cmap: str = "rainbow",
    figsize: tuple[float, float] | None = None,
    show: bool = True,
    save: str | Path | None = None,
    **kwargs: Any,
) -> tuple[plt.Figure, np.ndarray]:
    """Plot one gene across multiple spatial libraries using Squidpy."""
    apply_publication_defaults()
    import squidpy as sq

    if figsize is None:
        figsize = (5 * len(panels), 5)
    fig, axes = plt.subplots(1, len(panels), figsize=figsize, constrained_layout=True)
    axes_flat = as_axes_array(axes)

    for idx, lib in enumerate(panels):
        ax = sq.pl.spatial_scatter(
            adata,
            color=gene,
            layer=layer,
            library_key=library_key,
            library_id=[lib],
            title=[lib],
            spatial_key=spatial_key,
            size=size,
            img=False,
            frameon=False,
            cmap=cmap,
            ax=axes_flat[idx],
            return_ax=True,
            use_raw=False,
            **kwargs,
        )
        ax.set_aspect("equal")

    finish_figure(fig, save=save, show=show)
    return fig, axes_flat


def highlight_clusters_panels(
    adata: Any,
    color_key: str,
    select_cluster: Sequence[str],
    panels: Sequence[str],
    *,
    library_key: str = "library",
    spatial_key: str = "spatial",
    na_color: str = "#DEDEDE",
    size: float = 50,
    legend_title: str = "cluster",
    show: bool = True,
    save: str | Path | None = None,
    **kwargs: Any,
) -> tuple[plt.Figure, np.ndarray]:
    """Highlight selected spatial clusters across multiple libraries."""
    apply_publication_defaults()
    import squidpy as sq

    selected = {str(item) for item in select_cluster}
    adata_plt = adata.copy()
    cats = adata_plt.obs[color_key].astype("category").cat.categories
    color_map = pd.Series(adata_plt.uns[f"{color_key}_colors"], index=cats)
    adata_plt.uns[f"{color_key}_colors"] = [
        row if str(idx) in selected else na_color for idx, row in color_map.items()
    ]

    ncols = min(4, len(panels))
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5 * ncols, 5 * nrows),
        constrained_layout=True,
    )
    axes_flat = as_axes_array(axes)

    for idx, lib in enumerate(panels):
        ax = axes_flat[idx]
        sq.pl.spatial_scatter(
            adata_plt,
            color=color_key,
            library_key=library_key,
            library_id=[lib],
            title=[lib],
            spatial_key=spatial_key,
            size=size,
            img=False,
            na_color=na_color,
            frameon=False,
            legend_loc="on data",
            legend_fontsize=8,
            legend_fontweight="normal",
            ax=ax,
            return_ax=True,
            **kwargs,
        )
        ax.set_aspect("equal")

    for ax in axes_flat[len(panels):]:
        ax.axis("off")

    handles = [
        mpatches.Patch(color=color_map[cat], label=str(cat))
        for cat in sorted(cats, key=_try_int)
    ]
    fig.legend(
        handles=handles,
        loc="center right",
        bbox_to_anchor=(1.05, 0.5),
        ncol=1,
        frameon=False,
        title=legend_title,
    )
    finish_figure(fig, save=save, show=show)
    return fig, axes_flat


def highlight_and_expression_grid(
    adata: Any,
    row_panels: Sequence[str],
    column_panels: Sequence[str],
    group_key: str,
    *,
    select_cluster: Sequence[str] | None = None,
    select_cluster_highlight: bool = True,
    crop_to_selected: bool = False,
    library_key: str = "library",
    spatial_key: str = "spatial",
    spot_size: float = 50,
    crop_extend_factor: float = 5,
    na_color: str = "#DEDEDE",
    cmap: str = "rainbow",
    show: bool = True,
    save: str | Path | None = None,
    **kwargs: Any,
) -> tuple[plt.Figure, np.ndarray]:
    """Plot genes and cluster annotations as a spatial row-by-column grid."""
    apply_publication_defaults()
    import squidpy as sq

    adata_hl = adata.copy()
    cats = adata_hl.obs[group_key].astype("category").cat.categories
    color_map = pd.Series(adata_hl.uns[f"{group_key}_colors"], index=cats)
    selected = None if select_cluster is None else {str(item) for item in select_cluster}
    if selected is not None and select_cluster_highlight:
        adata_hl.uns[f"{group_key}_colors"] = [
            row if str(idx) in selected else na_color for idx, row in color_map.items()
        ]

    crop_coords: dict[str, list[tuple[int, int, int, int]]] = {}
    if crop_to_selected and selected is not None:
        extend = int(crop_extend_factor * spot_size)
        for lib in column_panels:
            mask = (
                adata.obs[group_key].astype(str).isin(selected)
                & (adata.obs[library_key] == lib)
            )
            coords = adata.obsm[spatial_key][mask.to_numpy()]
            if coords.shape[0] == 0:
                continue
            xmin, ymin = coords.min(axis=0)
            xmax, ymax = coords.max(axis=0)
            crop_coords[lib] = [
                (int(xmin) - extend, int(ymin) - extend, int(xmax) + extend, int(ymax) + extend)
            ]

    fig, axes = plt.subplots(
        len(row_panels),
        len(column_panels),
        figsize=(5 * len(column_panels), 5 * len(row_panels)),
        constrained_layout=True,
        squeeze=False,
    )

    for col_idx, lib in enumerate(column_panels):
        for row_idx, color_key in enumerate(row_panels):
            sq.pl.spatial_scatter(
                adata_hl,
                color=color_key,
                library_key=library_key,
                library_id=[lib],
                title=[lib],
                spatial_key=spatial_key,
                size=spot_size,
                crop_coord=crop_coords.get(lib),
                img=False,
                na_color=na_color,
                cmap=cmap,
                frameon=False,
                legend_loc="on data",
                legend_fontsize=8,
                legend_fontweight="normal",
                ax=axes[row_idx, col_idx],
                return_ax=True,
                **kwargs,
            )
            axes[row_idx, col_idx].set_aspect("equal")

    handles = [
        mpatches.Patch(color=color_map[cat], label=str(cat))
        for cat in sorted(cats, key=_try_int)
    ]
    fig.legend(
        handles=handles,
        loc="center right",
        bbox_to_anchor=(1.05, 0.5),
        ncol=1,
        frameon=False,
        title="cluster",
    )
    finish_figure(fig, save=save, show=show)
    return fig, axes
