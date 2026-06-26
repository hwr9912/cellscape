# 单细胞绘图

单细胞绘图建议通过子模块导入：

```python
import cellscape.single_cell as scp
```

## UMAP 类别与表达组合图

`umap_expr_with_category` 用于查看指定类别细胞中的基因表达，并在下一行展示该类别在
全局 UMAP 中的位置。

```python
fig = scp.umap_expr_with_category(
    adata,
    category_col="celltype",
    category_value=["alpha", "beta"],
    gene="GeneA",
    cmap="Reds",
    show=False,
)
```

常用参数：

- `category_col`：`adata.obs` 中的分类列。
- `category_value`：要展示的一个或多个类别；传入 `None` 时展示观测到的全部类别。
- `gene`：基因名，或 `adata.obs` 中的连续变量列名。
- `bg_color`：非高亮细胞的背景色。
- `legend_loc`：分类标签位置，默认 `"on data"`。

## 单个坐标轴高亮

如果你已经自己创建了 Matplotlib 画布，可以使用 `umap_highlight` 写入指定 `Axes`。

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(4, 4))

scp.umap_highlight(
    adata,
    ax,
    mode="category",
    category_col="celltype",
    category_value="alpha",
)
```

## Hotspot 局部相关性热图

`local_correlation_plot` 用于展示 Hotspot module 的局部相关性聚类结果。

```python
ax = scp.local_correlation_plot(
    local_correlation_z,
    modules,
    linkage,
    vmin=-8,
    vmax=8,
    show=False,
)
```

其中：

- `local_correlation_z`：局部相关性 Z-score 矩阵。
- `modules`：每个基因对应的 module 编号。
- `linkage`：层次聚类 linkage 结果。
