# cellscape

`cellscape` 是一个面向单细胞和空间组学分析的 Python 可视化与注释辅助工具包。
它以 `AnnData` 为主要输入，封装常见的空间表达面板、空间分类高亮、UMAP 高亮、
Hotspot 局部相关性热图，以及 labelme 区域注释到 mask / `adata.obs` 的转换流程。

## 适合做什么

- 按多个 library 对比空间表达量
- 在多个空间切片中高亮指定细胞类型或 cluster
- 把表达量和分类注释组合成 publication-ready 网格图
- 在 UMAP 上高亮某一类细胞，并查看该类中的基因表达
- 绘制 Hotspot module 局部相关性聚类热图
- 将 labelme 标注转换为 PNG / NPY mask，并投射回 `AnnData.obs`

## 安装

从 GitHub 安装：

```bash
pip install git+https://github.com/hwr9912/cellscape.git
```

本地开发安装：

```bash
pip install -e .
```

安装开发依赖：

```bash
pip install -e ".[dev]"
```

项目要求 Python `>=3.12`。主要依赖包括 `anndata`、`scanpy`、`squidpy`、
`matplotlib`、`numpy`、`pandas`、`scipy`、`seaborn`、`pillow`、`glasbey`
和 `tqdm`。

## 快速开始

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
    img=True,
    show=False,
)
fig.savefig("GeneA_panels.png", dpi=300, bbox_inches="tight")
```

## 空间组学绘图

绘制多个 library 中的基因或连续变量表达：

```python
fig, axes = spt.spatial_expression_panels(
    adata,
    gene="GeneA",
    panels=["library0", "library1", "library2"],
    library_key="library",
    layer="lognorm",
    size=30,
    cmap="rainbow",
    img=True,
    show=False,
)
```

高亮指定 cluster，并保留统一图例：

```python
fig, axes = spt.highlight_clusters_panels(
    adata,
    color_key="celltype",
    select_cluster=["Oligo"],
    panels=["library0", "library1", "library2"],
    library_key="library",
    img=True,
    show=False,
)
```

把多行表达量和最后一行分类注释组合成网格：

```python
fig, axes = spt.highlight_and_expression_grid(
    adata,
    row_panels=["GeneA", "GeneB"],
    column_panels=["library0", "library1"],
    group_key="celltype",
    select_cluster=["Oligo"],
    crop_to_selected=True,
    img_group=True,
    show=False,
)
```

## 单细胞绘图

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

如果已经在外部准备好 Hotspot 的局部相关性矩阵、module 标注和层次聚类结果，可以绘制
module 热图：

```python
ax = scp.local_correlation_plot(
    local_correlation_z,
    modules,
    linkage,
    show=False,
)
```

## labelme 注释流程

单个 labelme JSON 转 mask：

```python
import cellscape.datasets as dts

result = dts.labelme_to_mask(
    image_path="library0.png",
    annotation_json_path="library0.json",
    labels=["__ignore__", "_background_", "tumor"],
    output_dir="mask",
)

print(result.class_npy)
print(result.class_png_dir)
print(result.visualization_jpg)
```

批量转换：

```python
results = dts.labelme_to_masks(
    annotation_dir="annotations",
    image_dir="images",
    output_dir="mask",
    labels=["__ignore__", "_background_", "tumor"],
)
```

把 mask 投射回 `adata.obs`：

```python
adata = spt.project_labelme_masks_to_obs(
    adata,
    mask_path_dict={
        "library0": "mask/SegmentationClassNpy/library0.npy",
        "library1": "mask/SegmentationClassNpy/library1.npy",
    },
    labels=["__ignore__", "_background_", "tumor"],
    mask_format="npy",
    write_mode="merge",
    region_key="region",
    na_value="other",
    inplace=False,
)
```

## 文档托管到 GitHub Pages

本仓库已经包含 `mkdocs.yml` 和 `docs/` 目录，可以用 MkDocs 生成中文文档站点。

本地预览：

```bash
pip install mkdocs mkdocs-material
mkdocs serve
```

构建静态站点：

```bash
mkdocs build
```

部署到 GitHub Pages：

```bash
mkdocs gh-deploy
```

也可以在 GitHub 仓库的 **Settings -> Pages** 中配置 GitHub Actions 自动部署。

## 开发与测试

```bash
pip install -e ".[dev]"
pytest
ruff check .
```

绘图函数通常支持 `show` 和 `save` 参数。设置 `show=False` 适合在测试、批处理脚本或
notebook 中先拿到 `Figure` / `Axes`，再统一保存或展示。
