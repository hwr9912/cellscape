# cellscape 文档

`cellscape` 是用于单细胞和空间组学分析的 Python 绘图与注释辅助工具包。
它围绕 `AnnData` 数据结构，提供多 library 空间表达面板、分类高亮、UMAP 高亮、
Hotspot module 热图，以及 labelme 区域注释处理流程。

## 功能概览

- 空间表达面板：按多个 library 对比基因或连续变量表达。
- 空间分类高亮：在多个 library 中突出显示指定细胞类型或 cluster。
- 表达量与注释网格：将多行表达量和分类注释组合到同一张图中。
- 单细胞 UMAP 高亮：查看指定类别细胞的位置及其基因表达。
- Hotspot 热图：展示 module 局部相关性聚类结果。
- labelme 注释处理：将 JSON 标注转换成 mask，并投射回 `adata.obs`。

## 安装

```bash
pip install git+https://github.com/hwr9912/cellscape.git
```

本地开发：

```bash
pip install -e ".[dev]"
```

## 快速示例

```python
import scanpy as sc
import cellscape.spatial as spt

adata = sc.read_h5ad("example_spatial.h5ad")

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

继续阅读：

- [安装](installation.md)
- [快速开始](quickstart.md)
- [空间组学绘图](spatial.md)
- [单细胞绘图](single-cell.md)
- [labelme 注释流程](labelme.md)
- [GitHub Pages 部署](github-pages.md)
