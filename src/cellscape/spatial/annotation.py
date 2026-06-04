"""Spatial annotation projection helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Literal

import numpy as np
from PIL import Image
from tqdm.auto import tqdm

from cellscape.core.validation import get_spatial_coordinates


@contextmanager
def _pillow_image_size_limit(allow_large_images: bool) -> Iterator[None]:
    """临时控制 Pillow 的图片大小安全限制, 并在退出时恢复原值"""
    previous_limit = Image.MAX_IMAGE_PIXELS
    if allow_large_images:
        Image.MAX_IMAGE_PIXELS = None
    try:
        yield
    finally:
        Image.MAX_IMAGE_PIXELS = previous_limit


def _normalize_projection_labels(labels: Sequence[str]) -> tuple[str, ...]:
    """清理投射标签, 并拒绝空标签或重复标签"""
    normalized = tuple(str(label).strip() for label in labels if str(label).strip())
    if not normalized:
        raise ValueError("`labels` 至少需要包含一个有效标签")
    if len(set(normalized)) != len(normalized):
        raise ValueError("`labels` 中存在重复标签")
    return normalized


def _load_npy_masks(path: str | Path, labels: Sequence[str]) -> np.ndarray:
    """读取 labelme_to_masks 输出的 NPY mask, 并校验第一维与 labels 一致"""
    npy_path = Path(path)
    if not npy_path.exists():
        raise FileNotFoundError(f"找不到 npy mask 文件: {npy_path}")
    masks = np.asarray(np.load(npy_path))
    if masks.ndim != 3:
        raise ValueError(
            f"npy mask 必须是三维数组 `(n_labels, height, width)`, 当前形状为 {masks.shape}"
        )
    if masks.shape[0] != len(labels):
        raise ValueError(
            "npy 模式下 `labels` 长度必须和 npy 第一维长度相等: "
            f"labels={len(labels)}, npy 第一维={masks.shape[0]}"
        )
    return masks > 0


def _load_png_masks(path: str | Path, labels: Sequence[str]) -> np.ndarray:
    """读取 labelme_to_masks 输出的 PNG mask 目录, 并按文件名排序匹配 labels"""
    png_dir = Path(path)
    if not png_dir.exists():
        raise FileNotFoundError(f"找不到 png mask 目录: {png_dir}")
    if not png_dir.is_dir():
        raise ValueError(f"png 模式下每个 library 的路径必须是目录: {png_dir}")

    png_files = sorted(png_dir.glob("*.png"))
    if len(png_files) != len(labels):
        raise ValueError(
            "png 模式下 `labels` 长度必须和路径下 PNG 文件数相等: "
            f"labels={len(labels)}, PNG 文件数={len(png_files)}, 路径={png_dir}"
        )

    masks = []
    expected_shape: tuple[int, int] | None = None
    for png_file in png_files:
        with _pillow_image_size_limit(allow_large_images=True):
            with Image.open(png_file) as image:
                mask = np.asarray(image.convert("L")) > 0
        if expected_shape is None:
            expected_shape = mask.shape
        elif mask.shape != expected_shape:
            raise ValueError(
                f"同一 png mask 目录下所有文件尺寸必须一致: {png_dir}"
            )
        masks.append(mask)
    return np.stack(masks, axis=0)


def _validate_libraries(
    observed_libraries: Sequence[Any],
    mask_paths: Mapping[str, str | Path],
    *,
    library_key: str,
) -> tuple[str, ...]:
    """确认每个 adata library 都有对应的 mask 路径"""
    libraries = tuple(str(library) for library in dict.fromkeys(observed_libraries))
    missing = [library for library in libraries if library not in mask_paths]
    if missing:
        raise KeyError(
            f"`mask_paths` 缺少 adata.obs[{library_key!r}] 中的 library: {missing}"
        )
    return libraries


def _project_masks_to_rows(
    masks: np.ndarray,
    coords: np.ndarray,
) -> np.ndarray:
    """把坐标投射到多层 mask, 返回形状为 `(n_cells, n_labels)` 的布尔矩阵"""
    # 读取
    height, width = masks.shape[1], masks.shape[2]
    pixel_xy = np.rint(coords).astype(int)
    x = np.clip(pixel_xy[:, 0], 0, width - 1)
    y = np.clip(pixel_xy[:, 1], 0, height - 1)
    valid = (
        (coords[:, 0] >= 0)
        & (coords[:, 1] >= 0)
        & (coords[:, 0] <= width)
        & (coords[:, 1] <= height)
    )
    projected = np.zeros((coords.shape[0], masks.shape[0]), dtype=bool)
    if np.any(valid):
        # mask 数组按 [label, y, x] 存储, 空间坐标按 [x, y] 投射。
        projected[valid, :] = masks[:, y[valid], x[valid]].T
    return projected


def project_labelme_masks_to_obs(
    adata: Any,
    mask_path_dict: Mapping[str, str | Path],
    labels: Sequence[str],
    *,
    mask_format: Literal["npy", "png"] = "npy",
    write_mode: Literal["separate", "merge"] = "separate",
    region_key: str = "region",
    na_value: Any = "other",
    library_key: str = "library",
    spatial_key: str = "spatial",
    inplace: bool = False,
) -> Any | None:
    """
    将 labelme mask 图片注释按空间坐标投射到 AnnData 的 obs 布尔列

    该函数接受 `labelme_to_masks` 生成的 NPY 或 PNG mask。空间坐标按最近邻像素取整,
    默认写出为多列模式, 落在区域 mask 中的 cellbin 会在
    `adata.obs[f"in_{region}"]` 中标记为 True。也可以使用单列模式写入
    `adata.obs[region_key]`, 未命中的 cellbin 填入 `na_value`, 命中的 cellbin
    按 `labels` 顺序直接填入区域名; 如果多个区域重叠, 后写入的区域会覆盖先写入的区域。

    ### Parameters
    - `adata: Any` 包含 `obs[library_key]` 和 `obsm[spatial_key]` 的 AnnData 对象,
    - `mask_paths: Mapping[str, str | Path]` library 到 mask 路径的字典; NPY 模式为 npy 文件, PNG 模式为 PNG 目录,
    - `labels: Sequence[str]` 每层或每个 PNG 对应的区域名称,
    - `mask_format: Literal["npy", "png"] = "npy"` mask 输入格式,
    - `write_mode: Literal["separate", "merge"] = "separate"` 写出模式; `"separate"` 为每个区域写一个布尔列,
      `"merge"` 为所有区域写入 `obs[region_key]` 一列,
    - `region_key: str = "region"` 单列模式下写入 `adata.obs` 的列名,
    - `na_value: Any = "other"` 单列模式下未命中任一区域时填入的默认值,
    - `library_key: str = "library"` `adata.obs` 中记录 library 名称的列,
    - `spatial_key: str = "spatial"` `adata.obsm` 中记录空间坐标的键,
    - `inplace: bool = False` 是否直接修改输入 `adata`; 为 False 时返回修改后的 copy

    ### Example
    默认情况下可以直接使用
    ```python
    # 测试读取模式npy和写入模式separate
    import scanpy as sc
    import cellscape.spatial as spt

    adata = sc.read_h5ad("example_sp_data.h5ad")
    adata = spt.project_labelme_masks_to_obs(
        adata=adata,
        mask_path_dict={'lib0':"library_backgrounds/masks/SegmentationClassNpy/lib0.npy",
                        'lib1':"library_backgrounds/masks/SegmentationClassNpy/lib1.npy",
                        'lib2':"library_backgrounds/masks/SegmentationClassNpy/lib2.npy"},
        labels=["__ignore__", "_background_", "test"],
        inplace=False,
    )
    adata.obs["region"] = "other"
    adata.obs.loc[adata.obs["in_test"], "region"] = "test"
    ```
    也可以读取图片和合并输出到一列
    ```python
    # 测试读取模式png和写入模式merge
    import scanpy as sc
    import cellscape.spatial as spt

    adata = sc.read_h5ad("example_sp_data.h5ad")
    adata = spt.project_labelme_masks_to_obs(
        adata=adata,
        mask_path_dict={'lib0':"library_backgrounds/masks/SegmentationClass/lib0",
                        'lib1':"library_backgrounds/masks/SegmentationClass/lib1",
                        'lib2':"library_backgrounds/masks/SegmentationClass/lib2"},
        mask_format='png',
        labels=["__ignore__", "_background_", "test"],
        write_mode="merge",
        region_key="region",
        inplace=False,
    )
    ```
    """
    if mask_format not in {"npy", "png"}:
        raise ValueError("`mask_format` 只能是 'npy' 或 'png'")
    if write_mode not in {"separate", "merge"}:
        raise ValueError("`write_mode` 只能是 'separate' 或 'merge'")

    labels = _normalize_projection_labels(labels)
    if library_key not in adata.obs:
        raise KeyError(f"adata.obs 缺少 library 列: {library_key!r}")
    if spatial_key not in adata.obsm:
        raise KeyError(f"adata.obsm 缺少空间坐标键: {spatial_key!r}")
    coords = get_spatial_coordinates(adata, spatial_key=spatial_key)
    target = adata if inplace else adata.copy()

    libraries = _validate_libraries(
        target.obs[library_key].to_numpy(),
        mask_path_dict,
        library_key=library_key,
    )
    result = np.zeros((target.n_obs, len(labels)), dtype=bool)
    obs_libraries = target.obs[library_key].astype(str).to_numpy()

    for library in tqdm(libraries, desc="投射 library mask", unit="library"):
        row_mask = obs_libraries == library
        if mask_format == "npy":
            masks = _load_npy_masks(mask_path_dict[library], labels)
        else:
            masks = _load_png_masks(mask_path_dict[library], labels)
        result[row_mask, :] = _project_masks_to_rows(masks, coords[row_mask])

    if write_mode == "separate":
        for idx, region in tqdm(
            enumerate(labels),
            total=len(labels),
            desc="写入 obs 区域列",
            unit="region",
        ):
            target.obs.loc[:, f"in_{region}"] = result[:, idx]
    else:
        target.obs.loc[:, region_key] = na_value
        for idx, region in tqdm(
            enumerate(labels),
            total=len(labels),
            desc="写入 obs 区域列",
            unit="region",
        ):
            target.obs.loc[result[:, idx], region_key] = f"{region}"

    if inplace:
        return None
    return target
