# 空间组学绘图

空间绘图函数建议通过子模块导入：

```python
import cellscape.spatial as spt
```

## 多 library 表达面板

`spatial_expression_panels` 用于在多个 library 中绘制同一个基因或连续变量。

```python
fig, axes = spt.spatial_expression_panels(
    adata,
    gene="GeneA",
    panels=["library0", "library1", "library2"],
    library_key="library",
    spatial_key="spatial",
    layer="lognorm",
    size=30,
    cmap="rainbow",
    img=True,
    show=False,
)
```

常用参数：

- `gene`：基因名，或 `adata.obs` 中的连续变量列名。
- `panels`：要展示的 library 名称列表。
- `library_key`：`adata.obs` 中记录 library 的列名。
- `layer`：表达矩阵所在 layer；传入 `None` 时使用 `adata.X`。
- `img`：是否显示 squidpy 空间背景图。

## 多 library 分类高亮

`highlight_clusters_panels` 用于在多个 library 中高亮指定类别。

```python
fig, axes = spt.highlight_clusters_panels(
    adata,
    color_key="celltype",
    select_cluster=["Oligo", "Astro"],
    panels=["library0", "library1", "library2"],
    library_key="library",
    size=30,
    img=True,
    show=False,
)
```

如果 `select_cluster=None`，函数会展示全部类别的原始配色。

## 表达量与分类注释网格

`highlight_and_expression_grid` 会把 `row_panels` 中的基因或连续变量放在前几行，
并在最后一行追加 `group_key` 对应的分类注释。

```python
fig, axes = spt.highlight_and_expression_grid(
    adata,
    row_panels=["GeneA", "GeneB"],
    column_panels=["library0", "library1"],
    group_key="celltype",
    select_cluster=["Oligo"],
    select_cluster_highlight=True,
    crop_to_selected=True,
    crop_extend_factor=1,
    img_group=True,
    show=False,
)
```

常用场景：

- 比较多个基因在多个切片中的空间分布。
- 在最后一行保留细胞类型或 cluster 参照。
- 使用 `crop_to_selected=True` 聚焦到指定类别所在区域。
