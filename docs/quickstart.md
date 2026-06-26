# 快速开始

`cellscape` 的主要输入是 `AnnData`。空间绘图通常需要：

- `adata.obsm["spatial"]`：二维空间坐标。
- `adata.obs["library"]`：library 或切片 ID。
- `adata.obs["celltype"]`：细胞类型、cluster 或其他分类注释。
- `adata.layers["lognorm"]` 或 `adata.X`：表达矩阵。

## 空间表达面板

```python
import scanpy as sc
import cellscape.spatial as spt

adata = sc.read_h5ad("example_spatial.h5ad")

fig, axes = spt.spatial_expression_panels(
    adata,
    gene="GeneA",
    panels=["library0", "library1"],
    library_key="library",
    layer="lognorm",
    size=30,
    img=True,
    show=False,
)
fig.savefig("GeneA_panels.png", dpi=300, bbox_inches="tight")
```

## 空间分类高亮

```python
fig, axes = spt.highlight_clusters_panels(
    adata,
    color_key="celltype",
    select_cluster=["Oligo"],
    panels=["library0", "library1"],
    library_key="library",
    img=True,
    show=False,
)
```

## UMAP 表达与分类高亮

```python
import cellscape.single_cell as scp

fig = scp.umap_expr_with_category(
    adata,
    category_col="celltype",
    category_value=["alpha", "beta"],
    gene="GeneA",
    cmap="Reds",
    show=False,
)
```
