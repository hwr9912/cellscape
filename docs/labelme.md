# labelme 注释流程

`cellscape` 提供两个步骤来处理 labelme 区域注释：

1. 将 labelme JSON 转换为 PNG / NPY mask。
2. 按空间坐标把 mask 投射回 `adata.obs`。

## 单个 JSON 转 mask

```python
import cellscape.datasets as dts

result = dts.labelme_to_mask(
    image_path="images/library0.png",
    annotation_json_path="annotations/library0.json",
    labels=["__ignore__", "_background_", "tumor"],
    output_dir="mask",
)
```

返回值 `result` 记录了输出路径：

```python
print(result.class_npy)
print(result.class_png_dir)
print(result.visualization_jpg)
print(result.class_names_file)
```

## 批量转换

```python
results = dts.labelme_to_masks(
    annotation_dir="annotations",
    image_dir="images",
    output_dir="mask",
    labels=["__ignore__", "_background_", "tumor"],
)
```

输出目录通常包含：

```text
mask/
├── class_names.txt
├── SegmentationClass/
├── SegmentationClassNpy/
└── SegmentationClassVisualization/
```

## 投射到 AnnData

使用 NPY mask，并把每个区域写入独立布尔列：

```python
import cellscape.spatial as spt

adata = spt.project_labelme_masks_to_obs(
    adata,
    mask_path_dict={
        "library0": "mask/SegmentationClassNpy/library0.npy",
        "library1": "mask/SegmentationClassNpy/library1.npy",
    },
    labels=["__ignore__", "_background_", "tumor"],
    mask_format="npy",
    write_mode="separate",
    inplace=False,
)
```

这会写入：

```text
adata.obs["in___ignore__"]
adata.obs["in__background_"]
adata.obs["in_tumor"]
```

使用 PNG mask 目录，并合并写入单列：

```python
adata = spt.project_labelme_masks_to_obs(
    adata,
    mask_path_dict={
        "library0": "mask/SegmentationClass/library0",
        "library1": "mask/SegmentationClass/library1",
    },
    labels=["__ignore__", "_background_", "tumor"],
    mask_format="png",
    write_mode="merge",
    region_key="region",
    na_value="other",
    inplace=False,
)
```

## 注意事项

- `mask_path_dict` 的 key 必须能对应 `adata.obs[library_key]` 中的 library 名称。
- mask 图像尺寸应与空间坐标使用的原图坐标系一致。
- 如果多个区域重叠，`write_mode="merge"` 时后写入的区域会覆盖先写入的区域。
- 超大图像需要确认来源可信且内存足够后，再设置 `allow_large_images=True`。
