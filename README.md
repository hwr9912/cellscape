# cellscape

`cellscape` 是一个面向单细胞与空间组学分析的 Python 可视化工具包。当前版本以
`AnnData` 为主要输入，提供空间转录组绘图、单细胞 UMAP 高亮、Hotspot 局部相关性
热图，以及 labelme 区域注释到 mask / `adata.obs` 的转换工具。

## 主要功能

- 空间坐标散点图：支持 `adata.obsm["spatial"]`，也支持 `adata.obs` 中的 `x` / `y`
  坐标列。
- 空间表达面板：按多个 library 绘制基因或连续变量表达。
- 空间分类高亮：按 library 高亮指定 cluster，也可以保留全部分类原始配色。
- 表达量 + 分类注释网格：按行展示多个基因或连续变量，并在最后一行追加分类注释；
  支持按选中 cluster 裁剪视野。
- 细胞边界多边形绘图：读取 `adata.obsm["cell_border"]` 中的 cell polygon。
- 单细胞 UMAP 高亮：支持分类高亮、指定分类内基因表达，以及二者组合面板。
- Hotspot module 局部相关性聚类热图。
- labelme 注释预处理：把 labelme JSON 转成 PNG / NPY mask，并生成叠加预览图。
- labelme mask 投射：把区域 mask 按空间坐标写入 `adata.obs` 的布尔列或单个区域列。
- 统一的主题、调色板、校验和保存 / 展示收尾逻辑。

## 项目架构

```text
cellscape/
├── src/cellscape/
│   ├── __init__.py          # 顶层惰性导出，避免 import 时加载完整绘图依赖
│   ├── core/                # 共享绘图基础设施与 AnnData 输入校验
│   │   ├── plotting.py
│   │   └── validation.py
│   ├── spatial/             # 空间组学绘图与区域注释投射
│   │   ├── plots.py         # spatial_scatter / panel / grid 类绘图
│   │   ├── polygon.py       # cell_boundary_plot
│   │   └── annotation.py    # project_labelme_masks_to_obs
│   ├── single_cell/         # UMAP 与 Hotspot 绘图
│   │   └── plots.py
│   ├── datasets/            # labelme JSON 到 mask 的数据预处理工具
│   │   └── labelme.py
│   ├── styles/              # 主题和调色板
│   ├── io/                  # 读写辅助模块
│   └── widgets/             # 交互组件预留模块
├── examples/                # 最小用法示例
├── tests/                   # pytest 测试和 notebook 检查用例
├── docs/                    # 文档草稿与迁移记录
├── pyproject.toml           # hatchling 构建配置、依赖和工具配置
└── README.md
```

包采用 `src/` 布局，通过 `hatchling` 构建。顶层 `cellscape` API 使用
`__getattr__` 惰性导入绘图函数，因此 `import cellscape as cs` 不会立即执行分析、
读取数据或绘制图形。

## 安装

从 GitHub 安装：

```bash
pip install git+https://github.com/hwr9912/cellscape.git
```

本地开发安装：

```bash
pip install -e .
```

安装开发工具：

```bash
pip install -e ".[dev]"
```

可选的交互式依赖：

```bash
pip install -e ".[interactive]"
```

项目要求 Python `>=3.12`。核心依赖包括 `anndata`、`scanpy`、`squidpy`、
`matplotlib`、`numpy`、`pandas`、`scipy`、`seaborn`、`pillow`、`glasbey` 和
`tqdm`。

## 常用入口

可以从顶层包直接调用常用函数：

```python
import cellscape as cs

fig, ax = cs.spatial_scatter(adata, color="celltype", show=False)
fig, ax = cs.cell_boundary_plot(adata, color="cluster", show=False)
fig = cs.umap_expr_with_category(
    adata,
    category_col="celltype",
    category_value=["alpha", "beta"],
    gene="GeneA",
    show=False,
)
```

也可以按子模块导入，适合在脚本里明确区分工作流：

```python
import cellscape.spatial as spt
import cellscape.single_cell as scp
import cellscape.datasets as dts
```

### 空间组学绘图

```python
fig, ax = spt.spatial_scatter(
    adata,
    color="celltype",
    spatial_key="spatial",
    size=8,
    show=False,
)

fig, axes = spt.spatial_expression_panels(
    adata,
    gene="GeneA",
    panels=["library0", "library1"],
    library_key="library",
    layer="lognorm",
    img=True,
    show=False,
)

fig, axes = spt.highlight_clusters_panels(
    adata,
    color_key="celltype",
    select_cluster=["Oligo"],
    panels=["library0", "library1"],
    show=False,
)

fig, axes = spt.highlight_and_expression_grid(
    adata,
    row_panels=["GeneA", "GeneB"],
    column_panels=["library0", "library1"],
    group_key="celltype",
    select_cluster=["Oligo"],
    crop_to_selected=True,
    show=False,
)
```

### 细胞边界绘图

`cell_boundary_plot` 期望 `adata.obsm[border_key]` 的形状为
`(n_cells, n_vertices, 2)`。默认 `border_key="cell_border"`，并在
`relative=True` 时把 polygon 坐标平移到细胞中心坐标。

```python
fig, ax = spt.cell_boundary_plot(
    adata,
    color="celltype",
    border_key="cell_border",
    spatial_key="spatial",
    show=False,
)
```

### 单细胞绘图

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

fig = scp.umap_expr_with_category(
    adata,
    category_col="celltype",
    category_value=["alpha", "beta"],
    gene="GeneA",
    show=False,
)

ax = scp.local_correlation_plot(
    local_correlation_z,
    modules,
    linkage,
    show=False,
)
```

### labelme 注释和区域投射

单个 labelme JSON 转 mask：

```python
result = dts.labelme_to_mask(
    image_path="image.png",
    annotation_json_path="image.json",
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
    annotation_dir="library_backgrounds",
    image_dir="library_backgrounds",
    output_dir="library_backgrounds/mask",
    labels=["__ignore__", "_background_", "tumor"],
)
```

把 NPY mask 投射到 `adata.obs` 的多列布尔结果：

```python
adata = spt.project_labelme_masks_to_obs(
    adata,
    mask_path_dict={
        "library0": "library_backgrounds/mask/SegmentationClassNpy/library0.npy",
        "library1": "library_backgrounds/mask/SegmentationClassNpy/library1.npy",
    },
    labels=["__ignore__", "_background_", "tumor"],
    mask_format="npy",
    write_mode="separate",
    inplace=False,
)
```

把 PNG mask 目录合并投射到单个区域列：

```python
adata = spt.project_labelme_masks_to_obs(
    adata,
    mask_path_dict={
        "library0": "library_backgrounds/mask/SegmentationClass/library0",
        "library1": "library_backgrounds/mask/SegmentationClass/library1",
    },
    labels=["__ignore__", "_background_", "tumor"],
    mask_format="png",
    write_mode="merge",
    region_key="region",
    na_value="other",
    inplace=False,
)
```

## 开发与测试

```bash
pip install -e ".[dev]"
pytest
ruff check .
```

绘图函数在适用场景下支持 `show` 和 `save` 参数。设置 `show=False` 适合在测试、
批处理脚本或 notebook 中先拿到 `Figure` / `Axes`，再统一保存或展示。
