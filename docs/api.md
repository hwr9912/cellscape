# API 速查

## 空间组学

```python
import cellscape.spatial as spt
```

- `spt.spatial_expression_panels(...)`
- `spt.highlight_clusters_panels(...)`
- `spt.highlight_and_expression_grid(...)`
- `spt.project_labelme_masks_to_obs(...)`
- `spt.update_obs_from_df(...)`
- `spt.update_obs_from_bool_df(...)`

### obs 注释更新

```python
spt.update_obs_from_df(
    adata,
    df,
    index_columns,
    source_columns,
    target_columns,
    check_adata_index_unique=True,
    inplace=True,
)
```

`update_obs_from_df` 根据 `index_columns` 匹配 `df` 和 `adata.obs` 的行，并把
`df[source_columns]` 写入 `adata.obs[target_columns]`。`source_columns` 和
`target_columns` 必须同时是字符串，或同时是等长列表。

```python
spt.update_obs_from_bool_df(
    adata,
    df,
    index_columns,
    source_columns,
    target_columns,
    match_columns,
    match_values,
    update_values,
    check_adata_index_unique=True,
    inplace=True,
)
```

`update_obs_from_bool_df` 同样先按 `index_columns` 匹配行。只有当
`df[source_columns]` 为 True，且 `adata.obs[match_columns]` 当前值等于
`match_values` 时，才把 `adata.obs[target_columns]` 写为 `update_values`。

`source_columns`、`target_columns`、`match_columns`、`match_values` 和
`update_values` 必须同时是字符串，或同时是等长列表；`df[source_columns]`
必须是 bool 类型。

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
