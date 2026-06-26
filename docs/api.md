# API 速查

## 空间组学

```python
import cellscape.spatial as spt
```

- `spt.spatial_expression_panels(...)`
- `spt.highlight_clusters_panels(...)`
- `spt.highlight_and_expression_grid(...)`
- `spt.project_labelme_masks_to_obs(...)`

## 单细胞

```python
import cellscape.single_cell as scp
```

- `scp.umap_highlight(...)`
- `scp.umap_expr_with_category(...)`
- `scp.local_correlation_plot(...)`

## 数据集与注释预处理

```python
import cellscape.datasets as dts
```

- `dts.labelme_to_mask(...)`
- `dts.labelme_to_masks(...)`
- `dts.batch_labelme_to_masks(...)`
- `dts.LabelmeMaskResult`

## 统一参数

多数绘图函数支持：

- `show`：是否立即展示图像。
- `save`：保存路径。

在脚本或 notebook 中，推荐设置 `show=False`，拿到 `Figure` / `Axes` 后再统一保存。
