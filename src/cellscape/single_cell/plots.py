"""High-level plotting functions for single-cell analysis."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from cellscape.core.plotting import apply_publication_defaults, finish_figure
from cellscape.core.validation import require_obs_keys


def umap_highlight(
    adata: Any,
    ax: plt.Axes,
    *,
    mode: str = "category",
    category_col: str = "celltype",
    category_value: Any,
    gene: str | None = None,
    title: str | None = None,
    size: float = None,
    bg_color: str = "#d3d3d3",
    cmap: str | None = None,
    legend_loc: str = "on data",
    category_order: list[Any] | None = None,
    **kwargs: Any,
) -> plt.Axes:
    """显示所有细胞的同时高亮指定群
    - ``mode="category"`` 类别模式, 该模式下绘制类别图，高亮特定类别
    - ``mode="gene"`` 展示选定基因表达水平或者obs下某个连续变量的水平

    参数说明：
    - `adata: Any` 绘制依赖的数据,
    - `ax: plt.Axes` 传入作为画板的axes,
    - `mode: str = "category"` 绘制模式,
    - `category_col: str = "celltype"` 下面这个类别参数对应obs的那一列,
    - `category_value: Any` ,
        - 类别绘制模式下, 这个参数控制哪个类别高亮
        - 基因绘制模式下，这个参数控制显示哪个类别
    - `gene: str | None = None` (仅基因绘制模式下有效)基因或者obs中连续值,
    - `title: str | None = None` 标题, 对应sc.pl.umap中的同名参数,
    - `size: float = None` spot或cell的大小,
    - `bg_color: str = "#d3d3d3"` (仅类别绘制模式下有效)非高亮类别显示为什么颜色
    - `cmap: str | None = None` (仅基因绘制模式下有效)基因表达量图的连续颜色条,
    - `legend_loc: str = "on data"` (仅类别绘制模式下有效)图例显示位置,
    - `category_order: list[Any] | None = None` (仅类别绘制模式下有效, 一般用不到)颜色对应的类别顺序,
    """
    apply_publication_defaults()

    # 检查category_col参数是否在obs中
    require_obs_keys(adata, category_col)
    # 构建mask: 高亮 adata.obs[category_col] 中等于/不等于 category_value 的spot
    is_fg = (adata.obs[category_col] == category_value).to_numpy()
    is_bg = ~is_fg

    # # 如果没有任何 cell/spot 属于你想高亮的那个 group，就不画图，直接关闭坐标轴并返回 ax
    # # 这段作用存疑，先注释掉看一下
    # if not np.any(is_fg):
    #     ax.set_axis_off()
    #     return ax
    
    # ============ 1) 基因模式：只画当前组 ============
    if mode == "gene":
        if gene is None:
            raise ValueError("gene must be provided when mode='gene'.")
        sc.pl.umap(
            adata[is_fg],
            ax=ax,
            show=False,
            frameon=False,
            color=gene,
            cmap=cmap,
            title=title or str(gene),
            size=size,
            legend_loc="none",
            **kwargs,
        )
        return ax
    
    # ============ 2) 分类模式：背景层 + 高亮层 ============
    elif mode == "category":
        # 检查adata.obs[category_col]是否存在
        require_obs_keys(adata, category_col)

        # 生成类别排列顺序
        if category_order is None:
            if hasattr(adata.obs[category_col], "cat"):
                category_order = list(adata.obs[category_col].cat.categories)
            else:
                warnings.warn(
                    "类别顺序未指定且adata.obs[category_col].cat不存在, "
                    "按出现顺序生成列表"
                )
                category_order = list(
                    pd.Index(adata.obs[category_col].astype(str).unique())
                )
        else:
            category_order = list(category_order)

        # 确保对应的分类在adata里面有对应的颜色
        # 检查是否有类别顺序
        if hasattr(adata.obs[category_col], "cat"):
            cats = list(adata.obs[category_col].cat.categories)
        else:
            raise ValueError("adata.obs[category_col].cat.categories不存在, 无法对类别上色")
        if f"{category_col}_colors" in adata.uns:
            colors = list(adata.uns[f"{category_col}_colors"])
        else:
            raise ValueError("adata.uns没有对应类别的颜色")
        # 查找每个类别对应的颜色
        if len(colors) == len(cats):
            # 如果当前类别 cat 在颜色字典里，就取对应颜色；如果找不到，就用默认颜色：
            adata.uns[f"{category_col}_colors"] = [
                dict(zip(cats, colors)).get(cat, "#d3d3d3") for cat in category_order
            ]
        else:
            raise ValueError("类别数和颜色数不相等, 无法上色")
        
        # 绘制bg_color背景
        if np.any(is_bg):
            sc.pl.umap(
                adata[is_bg],
                ax=ax,
                show=False,
                frameon=False,
                color=None,
                title=title or str(category_value),
                size=size,
                legend_loc="none",
                **kwargs,
            )
            # 将绘制出来的整层散点对象全部改成 bg_color
            ax.collections[-1].set_color(bg_color)
            # 放置在z=0的高度上
            ax.collections[-1].set_zorder(0)
        # 绘制前景
        sc.pl.umap(
            adata[is_fg],
            ax=ax,
            show=False,
            frameon=False,
            color=category_col,
            title=title or str(category_value),
            size=size,
            legend_loc=legend_loc,
            **kwargs,
        )
        # 放置在z=1的高度上
        ax.collections[-1].set_zorder(1)

        return ax
    
    else:
        raise ValueError("参数 mode 必须是 'category' 或 'gene'.")


def umap_expr_with_category(
    adata: Any,
    *,
    category_col: str = "celltype",
    category_value: Any | list[Any] | None = None,
    gene: str | None = None,
    title_expr: str | None = None,
    title_cat: str | None = None,
    size: float = None,
    bg_color: str = "#d3d3d3",
    cmap: str | None = None,
    legend_loc: str = "on data",
    category_order: list[Any] | None = None,
    figsize: tuple[float, float] | None = None,
    show: bool = True,
    save: str | Path | None = None,
    **kwargs: Any,
) -> plt.Figure:
    """
    展示不同类型下的基因表达情况
    图片包括两行, ``len(category_value)`` 列，对于第 ``i`` 列有：
    - 第一行调用 ``umap_highlight`` 按 ``mode`` 展示当前类别的基因表达
      或分类结果
    - 第二行调用 ``umap_highlight`` 展示当前类别在所有细胞中的分类高亮位置

    示例:
    ```python
    import scanpy as sc
    import cellscape.single_cell as scp
    adata = sc.read_h5ad("example_sc_data.h5ad")
    scp.umap_expr_with_category(adata, 
                                category_col="celltype", 
                                category_value=["alpha", "beta"], 
                                gene="A",)
    ```

    参数说明:
    - `adata: Any` 绘制依赖的数据,
    - `category_col: str = "celltype"` 下面这个类别参数对应obs的那一列,
    - `category_value: Any | list[Any] | None = None` 显示哪些类别,
    - `gene: str | None = None` 基因或者obs中连续值,
    - `title_expr: str | None = None` 第一行基因表达图标题, 对应sc.pl.umap中的同名参数,
    - `title_cat: str | None = None` 第二行类别高亮图标题, 对应sc.pl.umap中的同名参数,
    - `size: float = None` spot或cell的大小,
    - `bg_color: str = "#d3d3d3"` 非高亮类别显示为什么颜色,
    - `cmap: str | None = None` 基因表达图的连续颜色条,
    - `legend_loc: str = "on data"` 图例显示位置,
    - `category_order: list[Any] | None = None` (一般用不到)颜色对应的类别顺序,
    - `figsize: tuple[float, float] | None = None`,
    - `show: bool = True`,
    - `save: str | Path | None = None`,
    """
    apply_publication_defaults()
    require_obs_keys(adata, category_col)

    # if mode not in {"category", "gene"}:
    #     raise ValueError("参数 mode 必须是 'category' 或 'gene'.")
    # if mode == "gene" and gene is None:
    #     raise ValueError("gene must be provided when mode='gene'.")

    if category_order is None and hasattr(adata.obs[category_col], "cat"):
        category_order = list(adata.obs[category_col].cat.categories)

    observed = pd.Index(adata.obs[category_col].dropna().unique())
    if category_value is None:
        if category_order is not None:
            category_value = [value for value in category_order if value in observed]
        elif hasattr(adata.obs[category_col], "cat"):
            category_value = list(observed)
        else:
            category_value = list(observed)
    elif isinstance(category_value, (str, bytes)):
        category_value = [category_value]
    else:
        category_value = list(category_value)

    if len(category_value) == 0:
        raise ValueError("category_value must contain at least one category.")

    missing = [value for value in category_value if value not in observed]
    if missing:
        raise ValueError(
            "category_value contains values not present in "
            f"adata.obs[{category_col!r}]: "
            f"{missing}"
        )

    if figsize is None:
        figsize = (5 * len(category_value), 10)

    fig, axes = plt.subplots(
        2,
        len(category_value),
        figsize=figsize,
        constrained_layout=True,
        squeeze=False,
    )

    for idx, value in enumerate(category_value):
        # 基因表达
        umap_highlight(
            adata,
            axes[0, idx],
            mode="gene",
            category_col=category_col,
            category_value=value,
            gene=gene,
            title=title_expr or f"{value}: {gene}",
            size=size,
            bg_color=bg_color,
            cmap=cmap,
            legend_loc=legend_loc,
            category_order=category_order,
            **kwargs,
        )
        # 类别高亮
        umap_highlight(
            adata,
            axes[1, idx],
            mode="category",
            category_col=category_col,
            category_value=value,
            title=title_cat,
            size=size,
            bg_color=bg_color,
            legend_loc=legend_loc,
            category_order=category_order,
            **kwargs,
        )

    finish_figure(fig, save=save, show=show)
    return fig


def local_correlation_plot(
    local_correlation_z: pd.DataFrame,
    modules: pd.Series,
    linkage: Any,
    *,
    mod_cmap: str = "tab20",
    vmin: float = -8,
    vmax: float = 8,
    z_cmap: str = "RdBu_r",
    yticklabels: bool = False,
    show: bool = False,
    save: str | Path | None = None,
) -> plt.Axes:
    """Plot a Hotspot module local-correlation heatmap."""
    apply_publication_defaults()
    import seaborn as sns
    from scipy.cluster.hierarchy import leaves_list

    colors = list(plt.get_cmap(mod_cmap).colors)
    module_colors = {i: colors[(i - 1) % len(colors)] for i in modules.unique()}
    module_colors[-1] = "#ffffff"

    row_colors = pd.DataFrame(
        {"Modules": pd.Series([module_colors[i] for i in modules], index=local_correlation_z.index)}
    )
    cm = sns.clustermap(
        local_correlation_z,
        row_linkage=linkage,
        col_linkage=linkage,
        vmin=vmin,
        vmax=vmax,
        cmap=z_cmap,
        xticklabels=False,
        yticklabels=yticklabels,
        row_colors=row_colors,
        rasterized=True,
    )

    ax = cm.ax_heatmap
    ax.set_ylabel("")
    ax.set_xlabel("")
    cm.ax_row_dendrogram.remove()

    order = leaves_list(linkage)
    mod_reordered = modules.iloc[order]
    y_positions = np.arange(modules.size)
    for mod in mod_reordered.unique():
        if mod == -1:
            continue
        cm.ax_row_colors.text(
            -0.5,
            y=y_positions[mod_reordered == mod].mean(),
            s=f"Module {mod}",
            horizontalalignment="right",
            verticalalignment="center",
        )
    cm.ax_row_colors.set_xticks([])

    if cm.cax is not None:
        cm.cax.set_ylabel("Z-Scores")
        cm.cax.yaxis.set_label_position("left")

    finish_figure(cm.fig, save=save, show=show)
    return ax
