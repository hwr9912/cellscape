"""High-level plotting functions for spatial omics analysis."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any
import warnings

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import squidpy as sq

from cellscape.core.plotting import (
    apply_publication_defaults,
    as_axes_array,
    finish_figure,
)
from cellscape.core.validation import get_spatial_coordinates, require_color_key
from cellscape.styles.palettes import glasbey_palette


def _as_list(value: str | Sequence[str]) -> list[str]:
    if isinstance(value, str):
        return [value]
    return list(value)


def _try_int(value: Any) -> tuple[int, int | str]:
    try:
        return (0, int(value))
    except (TypeError, ValueError):
        return (1, str(value))


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
    *, # 它后面的参数必须用关键字传参，不能再用位置传参。
    library_key: str = "library",
    spatial_key: str = "spatial",
    layer: str | None = "lognorm",
    size: float = 50,
    cmap: str = "rainbow",
    figsize: tuple[float, float] | None = None,
    img: str | bool = False,
    show: bool = True,
    save: str | Path | None = None,
    **kwargs: Any,
) -> tuple[plt.Figure, np.ndarray]:
    """
    在选定的多个library下使用squidpy绘制某个基因的表达量

    ⚠️注意: 该函数要求必须有多个library, 否则会报错

    ### Parameters
    - `adata: Any` 绘制依赖的数据,
    - `gene: str` 基因或者obs中连续值列,
    - `panels: Sequence[str]` 指定library的名称列表,
    - `library_key: str = "library"` 存储library的名称对应的obs列名,
    - `spatial_key: str = "spatial"` 对应`sq.pl.spatial_scatter`同名参数,
    - `layer: str | None = "lognorm"` 使用anndata的哪个layer,
    - `size: float = 50` spot大小,
    - `cmap: str = "rainbow"` 连续颜色条,
    - `figsize: tuple[float, float] | None = None` 整个fig的尺寸, 对于None默认每个小图是5*5的正方形,
    - `img: str | bool = False` 是否显示背景图,
    - `show: bool = True` 是否显示,
    - `save: str | Path | None = None` 保存位置

    ### Example
    ```python
    import scanpy as sc
    import matplotlib.pyplot as plt
    import cellscape.spatial as spt

    adata = sc.read_h5ad("example_sp_data.h5ad")
    fig, axes = spt.spatial_expression_panels(
        adata=adata,
        gene="Gene10",
        panels=['lib1','lib0'],
        library_key='library',
        spatial_key='spatial',
        layer='lognorm',
        size=20,
        img=True
    )
    plt.show()
    ```
    """
    apply_publication_defaults()

    if figsize is None:
        figsize = (5 * len(panels), 5)
    fig, axes = plt.subplots(1, len(panels), figsize=figsize, constrained_layout=True)
    axes_flat = axes.flatten()

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
            img=img,
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
    select_cluster: str | Sequence[str] | None,
    panels: Sequence[str],
    *,
    library_key: str = "library",
    spatial_key: str = "spatial",
    na_color: str = "#DEDEDE",
    size: float = 50,
    legend_title: str = "cluster",
    ncols: int | None = None,
    img: str | bool = False,
    show: bool = True,
    save: str | Path | None = None,
    **kwargs: Any,
) -> tuple[plt.Figure, np.ndarray]:
    """
    对于多library数据高亮指定的群

    ### Parameters:
    - `adata: Any` 绘制依赖的数据,
    - `color_key: str` 选择哪个obs列作为上色依据,
    - `select_cluster: str | Sequence[str] | None` 选择哪些/哪个值作为高亮群, 为None时显示所有群的原始配色,
    - `panels: Sequence[str]` 每个panel显示哪个library,
    - `library_key: str = "library"` 存储library的名称对应的obs列名,
    - `spatial_key: str = "spatial"` 对应`sq.pl.spatial_scatter`同名参数,
    - `na_color: str = "#DEDEDE"` 指定上色obs列有NA值对应的颜色,
    - `size: float = 50` spot的大小,
    - `legend_title: str = "cluster"` 图例的标题,
    - `ncols: int | None = None` panel列数, 为None时最多每行4列,
    - `img: str | bool = False` 是否显示背景图,
    - `show: bool = True` 是否显示,
    - `save: str | Path | None = None` 保存位置

    ### Example:
    ```python
    import scanpy as sc
    import matplotlib.pyplot as plt
    import cellscape.spatial as spt

    adata = sc.read_h5ad("example_sp_data.h5ad")
    fig, axes = spt.highlight_clusters_panels(
        adata=adata,
        color_key="celltype",
        select_cluster="Oligo",
        panels=['lib0','lib1','lib2'],
        library_key='library',
        spatial_key='spatial',
        layer='lognorm',
        size=20,
        img=True,
        show=True,
    )
    ```
    """
    apply_publication_defaults()

    # 提取所有类别
    adata_plt = adata.copy()
    cats = adata_plt.obs[color_key].astype("category").cat.categories
    # 提取颜色，没有预设颜色时按类别数量自动生成
    color_key_colors = f"{color_key}_colors"
    if (
        color_key_colors not in adata_plt.uns
        or len(adata_plt.uns[color_key_colors]) != len(cats)
    ):
        warnings.warn(
            f"adata.uns['{color_key_colors}'] is missing or has a length "
            f"different from {color_key!r} categories; generating {len(cats)} "
            "glasbey colors.",
            stacklevel=2,
        )
        adata_plt.uns[color_key_colors] = glasbey_palette(len(cats))
    color_map = pd.Series(adata_plt.uns[color_key_colors], index=cats)
    if select_cluster is not None:
        # 提取高亮群
        selected = {str(item) for item in _as_list(select_cluster)}
        adata_plt.uns[color_key_colors] = [
            row if str(idx) in selected else na_color for idx, row in color_map.items()
        ]

    if ncols is None:
        ncols = min(4, len(panels))
    if ncols < 1:
        raise ValueError("ncols must be at least 1.")
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
            img=img,
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
    
    # 合并图例
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
    # 保存图片
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
    crop_extend_factor: float = 1,
    na_color: str = "#DEDEDE",
    cmap: str = "rainbow",
    img_group: str | bool = False,
    show: bool = True,
    save: str | Path | None = None,
    **kwargs: Any,
) -> tuple[plt.Figure, np.ndarray]:
    """
    绘制多library下的表达量和分组注释网格图

    该函数会把 `row_panels` 中的基因或obs连续变量作为前几行，
    并在最后额外追加一行 `group_key` 分类注释图。

    ### Parameters:
    - `adata: Any` 绘制依赖的数据,
    - `row_panels: Sequence[str]` 每一行显示的基因或obs列名,
    - `column_panels: Sequence[str]` 每一列显示的library名称,
    - `group_key: str` 用于最后一行分类注释和高亮的obs列名,
    - `select_cluster: Sequence[str] | None = None` 需要高亮或裁剪的分类值,
    - `select_cluster_highlight: bool = True` 是否只保留选中分类的原颜色，其他分类置为 `na_color`,
    - `crop_to_selected: bool = False` 是否根据选中分类的空间坐标裁剪每个library,
    - `library_key: str = "library"` 存储library名称的obs列名,
    - `spatial_key: str = "spatial"` 对应`sq.pl.spatial_scatter`同名参数,
    - `spot_size: float = 50` spot大小,
    - `crop_extend_factor: float = 1` 裁剪边界外扩系数，实际外扩距离为 `crop_extend_factor * spot_size`,
        - 这一项注意不要设置的过大, 否则会让图片失焦
    - `na_color: str = "#DEDEDE"` 未选中分类或NA值对应的颜色,
    - `cmap: str = "rainbow"` 连续变量使用的颜色条,
    - `img_group: str | bool = False` 最后一行 `group_key` 分类图是否显示背景图,
    - `show: bool = True` 是否显示,
    - `save: str | Path | None = None` 保存位置,
    - `**kwargs: Any` 传递给`sq.pl.spatial_scatter`的其他参数

    ### Example:
    ```python
    import scanpy as sc
    import matplotlib.pyplot as plt
    import cellscape.spatial as spt

    adata = sc.read_h5ad("example_sp_data.h5ad")
    fig, axes = spt.highlight_and_expression_grid(
        adata=adata,
        row_panels=["Gene10", "Gene14"],
        column_panels=["lib0", "lib1", "lib2"],
        group_key="celltype",
        select_cluster="Oligo",
        select_cluster_highlight=True,
        crop_to_selected=True,
        crop_extend_factor=1,
        img_group=True,
        show=True,
    )
    plt.show()
    ```
    """
    # matplotlib全局绘图常数
    apply_publication_defaults()
    # 硬复制防篡改原始数据
    adata_plt = adata.copy()
    # 提取最下面一行分组变量的所有类别
    cats = list(adata_plt.obs[group_key].astype("category").cat.categories)
    # 提取颜色，没有预设颜色时按类别数量自动生成
    group_key_colors = f"{group_key}_colors"
    if (
        group_key_colors not in adata_plt.uns
        or len(adata_plt.uns[group_key_colors]) != len(cats)
    ):
        warnings.warn(
            f"adata.uns['{group_key_colors}'] is missing or has a length "
            f"different from {group_key!r} categories; generating {len(cats)} "
            "glasbey colors.",
            stacklevel=2,
        )
        adata_plt.uns[group_key_colors] = glasbey_palette(len(cats))
    color_map = pd.Series(adata_plt.uns[group_key_colors], index=cats)
    # 只保留指定 cluster 的原颜色，其他 cluster 全部改成 na_color
    selected = None if select_cluster is None else {str(item) for item in _as_list(select_cluster)}
    if selected is not None and select_cluster_highlight:
        adata_plt.uns[group_key_colors] = [
            row if str(idx) in selected else na_color for idx, row in color_map.items()
        ]
    # 针对指定cluster裁剪坐标存储为crop_coords
    crop_coords: dict[str, list[tuple[int, int, int, int]]] = {}
    if crop_to_selected and selected is not None:
        # 根据系数计算外延范围
        extend = int(crop_extend_factor * spot_size)
        # 逐library遍历
        for lib in column_panels:
            # 计算上下界
            mask = (
                adata.obs[group_key].astype(str).isin(selected)
                & (adata.obs[library_key] == lib)
            )
            coords = adata.obsm[spatial_key][mask.to_numpy()]
            if coords.shape[0] == 0:
                continue
            xmin, ymin = coords.min(axis=0)
            xmax, ymax = coords.max(axis=0)
            # 基于外延范围更新上下界
            crop_coords[lib] = [
                (int(xmin) - extend, int(ymin) - extend, int(xmax) + extend, int(ymax) + extend)
            ]

    # 制图部分
    # 最下面加一行group_key
    fig, axes = plt.subplots(
        len(row_panels) + 1,
        len(column_panels),
        figsize=(5 * len(column_panels), 5 * (len(row_panels) + 1)),
        constrained_layout=True,
        squeeze=False,
    )

    all_row_panels = list(row_panels) + [group_key]
    for col_idx, lib in enumerate(column_panels):
        for row_idx, color_key in enumerate(all_row_panels):
            sq.pl.spatial_scatter(
                adata_plt,
                color=color_key,
                library_key=library_key,
                library_id=[lib],
                title=[lib],
                spatial_key=spatial_key,
                size=spot_size,
                crop_coord=crop_coords.get(lib),
                img=img_group if color_key == group_key else False,
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
        loc="lower right",
        bbox_to_anchor=(1.05, 0.02),
        ncol=1,
        frameon=False,
        title="cluster",
    )
    finish_figure(fig, save=save, show=show)
    return fig, axes
