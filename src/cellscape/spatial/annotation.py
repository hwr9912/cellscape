"""Spatial annotation projection helpers."""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Literal

import numpy as np
import pandas as pd
from PIL import Image
from tqdm.auto import tqdm

from cellscape.core.validation import get_spatial_coordinates


def _as_column_tuple(columns: str | Sequence[str], *, name: str) -> tuple[str, ...]:
    """把单个列名或列名序列标准化为 tuple, 并拒绝空值。"""
    if isinstance(columns, str):
        normalized = (columns,)
    else:
        normalized = tuple(columns)
    if not normalized:
        raise ValueError(f"`{name}` 至少需要包含一个列名")
    if not all(isinstance(column, str) and column for column in normalized):
        raise ValueError(f"`{name}` 必须是非空字符串或非空字符串列表")
    return normalized


def _column_arg_kind(value: str | Sequence[str], *, name: str) -> str:
    """返回列参数是单个字符串还是列表式参数。"""
    if isinstance(value, str):
        return "str"
    if isinstance(value, Sequence):
        return "sequence"
    raise ValueError(f"`{name}` 必须是非空字符串或非空字符串列表")


def _normalize_paired_columns(
    source_columns: str | Sequence[str],
    target_columns: str | Sequence[str],
) -> tuple[tuple[str, ...], tuple[str, ...], str]:
    """标准化来源列和目标列, 并要求二者形态一致。"""
    source_kind = _column_arg_kind(source_columns, name="source_columns")
    target_kind = _column_arg_kind(target_columns, name="target_columns")
    if source_kind != target_kind:
        raise ValueError(
            "`source_columns` 和 `target_columns` 必须同时是字符串, "
            "或同时是等长字符串列表"
        )

    normalized_source_columns = _as_column_tuple(
        source_columns,
        name="source_columns",
    )
    normalized_target_columns = _as_column_tuple(
        target_columns,
        name="target_columns",
    )
    if len(normalized_source_columns) != len(normalized_target_columns):
        raise ValueError(
            "`source_columns` 和 `target_columns` 必须同时是字符串, "
            "或同时是等长字符串列表"
        )
    return normalized_source_columns, normalized_target_columns, source_kind


def _normalize_bool_mode_values(
    match_values: str | Sequence[str] | None,
    update_values: str | Sequence[str] | None,
    *,
    expected_kind: str,
    expected_length: int,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """标准化 bool 更新的匹配值和写入值。"""
    if match_values is None:
        raise ValueError("`match_values` 是 bool 更新的必需参数")
    if update_values is None:
        raise ValueError("`update_values` 是 bool 更新的必需参数")

    match_kind = _column_arg_kind(match_values, name="match_values")
    update_kind = _column_arg_kind(update_values, name="update_values")
    if match_kind != expected_kind or update_kind != expected_kind:
        raise ValueError(
            "bool 更新中 `source_columns`, `target_columns`, `match_values` "
            "和 `update_values` 必须同时是字符串, 或同时是等长列表"
        )

    normalized_match_values = _as_column_tuple(match_values, name="match_values")
    normalized_update_values = _as_column_tuple(update_values, name="update_values")
    if (
        len(normalized_match_values) != expected_length
        or len(normalized_update_values) != expected_length
    ):
        raise ValueError(
            "bool 更新中 `source_columns`, `target_columns`, `match_values` "
            "和 `update_values` 必须同时是字符串, 或同时是等长列表"
        )
    return normalized_match_values, normalized_update_values


def _duplicated_key_preview(frame: pd.DataFrame, columns: Sequence[str]) -> list[Any]:
    """返回重复索引键的少量预览, 用于错误信息。"""
    duplicated = frame.loc[:, list(columns)].duplicated(keep=False)
    preview = frame.loc[duplicated, list(columns)].drop_duplicates().head(5)
    if len(columns) == 1:
        return preview.iloc[:, 0].tolist()
    return [tuple(row) for row in preview.to_numpy()]


def _ensure_unique_keys(
    frame: pd.DataFrame,
    columns: Sequence[str],
    *,
    frame_name: str,
) -> None:
    """检查指定列组合是否唯一。"""
    if frame.loc[:, list(columns)].duplicated(keep=False).any():
        preview = _duplicated_key_preview(frame, columns)
        raise ValueError(
            f"{frame_name} 中索引列 {list(columns)!r} 对应的索引不唯一; "
            f"重复示例: {preview}"
        )


def _ensure_assignable_categories(
    obs: pd.DataFrame,
    column: str,
    values: pd.Series,
) -> None:
    """如果目标列是 category, 先补齐即将写入的新类别。"""
    if not isinstance(obs[column].dtype, pd.CategoricalDtype):
        return
    new_values = pd.Series(values).dropna().unique()
    categories = obs[column].cat.categories
    missing = [value for value in new_values if value not in categories]
    if missing:
        obs[column] = obs[column].cat.add_categories(missing)


def update_obs_from_df(
    adata: Any,
    df: pd.DataFrame,
    index_columns: str | Sequence[str],
    source_columns: str | Sequence[str],
    target_columns: str | Sequence[str],
    *,
    check_adata_index_unique: bool = True,
    inplace: bool = True,
) -> Any | None:
    """
    使用外部 dataframe 按索引列局部更新 `adata.obs` 中的列。

    函数会根据 `index_columns` 在 `df` 和 `adata.obs` 之间匹配行。
    能在 `df` 中找到对应索引的 `adata.obs` 行会被更新;
    找不到对应索引的行保持不变。
    `source_columns` 是 `df` 中提供新值的列, `target_columns` 是
    `adata.obs` 中被写入的目标列, 二者必须同时是字符串或同时是等长列表。

    ### Parameters
    - `adata: Any` 原始 AnnData 对象
    - `df: pd.DataFrame` 包含索引列和替换列的 dataframe
    - `index_columns: str | Sequence[str]` 用于匹配的索引列名
    - `source_columns: str | Sequence[str]` `df` 中提供新值的列名
    - `target_columns: str | Sequence[str]` `adata.obs` 中被更新的列名
    - `check_adata_index_unique: bool = True` 是否把 `adata.obs`
      索引不唯一作为错误; 为 False 时仍检查, 但只发出 warning
    - `inplace: bool = True` 是否原位修改输入 `adata`
    """
    # 统一把单个索引列名和多个索引列名都转换成 tuple, 后续逻辑就可以
    # 用同一种方式处理单列键和多列组合键。
    index_columns = _as_column_tuple(index_columns, name="index_columns")

    # 标准化来源列和目标列, 并要求二者同形态等长, 确保每个 df 来源列
    # 都能精确对应一个 adata.obs 目标列。
    source_columns, target_columns, _ = _normalize_paired_columns(
        source_columns,
        target_columns,
    )

    # df 必须同时包含用于匹配行的索引列, 以及提供新值的来源列。
    # 索引列缺失时无法定位行; 来源列缺失时无法生成写入值。
    missing_df_index_columns = [column for column in index_columns if column not in df]
    if missing_df_index_columns:
        raise KeyError(f"df 缺少索引列: {missing_df_index_columns}")

    missing_source_columns = [column for column in source_columns if column not in df]
    if missing_source_columns:
        raise KeyError(f"df 缺少来源列: {missing_source_columns}")

    # adata.obs 也必须有同名索引列, 否则无法和 df 建立行级对应关系。
    missing_obs_index_columns = [
        column for column in index_columns if column not in adata.obs
    ]
    if missing_obs_index_columns:
        raise KeyError(f"adata.obs 缺少索引列: {missing_obs_index_columns}")

    # df 作为更新来源, 每个索引键只能出现一次; 否则同一个 obs 行会对应多条
    # 来源记录, 无法确定应该采用哪一个值。
    _ensure_unique_keys(df, index_columns, frame_name="df")

    # adata.obs 中重复索引键的处理由 check_adata_index_unique 控制。严格模式
    # 直接报错; 宽松模式允许一条 df 记录更新多个 obs 行, 但发出 warning。
    duplicated_obs_keys = adata.obs.loc[:, list(index_columns)].duplicated(keep=False)
    if duplicated_obs_keys.any():
        preview = _duplicated_key_preview(adata.obs, index_columns)
        message = (
            f"adata.obs 中索引列 {list(index_columns)!r} 对应的索引不唯一; "
            f"重复示例: {preview}"
        )
        if check_adata_index_unique:
            raise ValueError(message)
        warnings.warn(message, UserWarning, stacklevel=2)

    # inplace=True 时直接修改传入对象; inplace=False 时先复制, 最后返回副本。
    target = adata if inplace else adata.copy()

    # 如果目标列不存在, 先创建 NA 列。这样匹配行可以写入 df 的新值,
    # 未匹配行保持 NA。
    missing_target_columns = [
        column for column in target_columns if column not in target.obs
    ]
    for column in missing_target_columns:
        target.obs.loc[:, column] = pd.NA
    if missing_target_columns:
        warnings.warn(
            "adata.obs 缺少目标列, 已创建并填充为 NA: "
            f"{missing_target_columns}",
            UserWarning,
            stacklevel=2,
        )

    # 用 MultiIndex 表示单列或多列组合键, 这样单索引和复合索引可以共享同一套
    # 匹配逻辑。matched 是一个与 adata.obs 等长的布尔数组。
    df_keys = pd.MultiIndex.from_frame(df.loc[:, list(index_columns)])
    obs_keys = pd.MultiIndex.from_frame(target.obs.loc[:, list(index_columns)])
    matched = obs_keys.isin(df_keys)
    if not matched.any():
        raise ValueError(
            f"adata.obs 中没有任何行能根据索引列 {list(index_columns)!r} 匹配到 df"
        )

    # 逐对处理来源列和目标列: 对所有匹配行, 直接用 df 来源列的值覆盖
    # adata.obs 的目标列。
    for source_column, target_column in zip(
        source_columns,
        target_columns,
        strict=True,
    ):
        # 把 df 的来源列转成以 df_keys 为索引的 Series, 后续通过 reindex 按
        # obs_keys 顺序取值, 保证写入值和 adata.obs 的匹配行顺序一致。
        values_by_key = pd.Series(df[source_column].to_numpy(), index=df_keys)
        new_values = values_by_key.reindex(obs_keys[matched])

        # 如果目标列是 category, 先补齐新类别, 避免 pandas 赋值时报错。
        _ensure_assignable_categories(target.obs, target_column, new_values)
        target.obs.loc[matched, target_column] = new_values.to_numpy()

    # 原位模式遵循 pandas/anndata 常见习惯返回 None; 非原位模式返回更新副本。
    if inplace:
        return None
    return target


def update_obs_from_bool_df(
    adata: Any,
    df: pd.DataFrame,
    index_columns: str | Sequence[str],
    source_columns: str | Sequence[str],
    target_columns: str | Sequence[str],
    match_values: str | Sequence[str],
    update_values: str | Sequence[str],
    *,
    check_adata_index_unique: bool = True,
    inplace: bool = True,
) -> Any | None:
    """
    使用 dataframe 中的 bool 列按条件更新 `adata.obs` 中的列。

    函数会先根据 `index_columns` 在 `df` 和 `adata.obs` 之间匹配行。
    对匹配行, 只有当 `df[source_columns]` 为 True, 且
    `adata.obs[target_columns]` 当前值等于 `match_values` 时, 才把目标列
    写为 `update_values`。

    `source_columns`, `target_columns`, `match_values`, `update_values` 必须
    同时是字符串, 或同时是等长字符串列表。

    ### Parameters
    - `adata: Any` 原始 AnnData 对象
    - `df: pd.DataFrame` 包含索引列和 bool 来源列的 dataframe
    - `index_columns: str | Sequence[str]` 用于匹配的索引列名
    - `source_columns: str | Sequence[str]` `df` 中作为更新开关的 bool 列名
    - `target_columns: str | Sequence[str]` `adata.obs` 中被更新的列名
    - `match_values: str | Sequence[str]` 目标列当前需要匹配的值
    - `update_values: str | Sequence[str]` 条件满足时写入目标列的值
    - `check_adata_index_unique: bool = True` 是否把 `adata.obs`
      索引不唯一作为错误; 为 False 时仍检查, 但只发出 warning
    - `inplace: bool = True` 是否原位修改输入 `adata`
    """
    # 统一索引列形态, 让单列键和多列组合键共用后续匹配逻辑。
    index_columns = _as_column_tuple(index_columns, name="index_columns")

    # bool 更新涉及四组一一对应的参数: df bool 来源列、obs 目标列、
    # 目标列当前匹配值、条件满足后的写入值。这里先校验列参数同形态等长。
    source_columns, target_columns, column_kind = _normalize_paired_columns(
        source_columns,
        target_columns,
    )
    match_values, update_values = _normalize_bool_mode_values(
        match_values,
        update_values,
        expected_kind=column_kind,
        expected_length=len(source_columns),
    )

    # df 必须包含索引列和 bool 来源列; 来源列必须是真正的 bool dtype,
    # 因为这些列会被解释为“是否允许更新”的开关。
    missing_df_index_columns = [column for column in index_columns if column not in df]
    if missing_df_index_columns:
        raise KeyError(f"df 缺少索引列: {missing_df_index_columns}")

    missing_source_columns = [column for column in source_columns if column not in df]
    if missing_source_columns:
        raise KeyError(f"df 缺少来源列: {missing_source_columns}")

    non_bool_source_columns = [
        column
        for column in source_columns
        if not pd.api.types.is_bool_dtype(df[column].dtype)
    ]
    if non_bool_source_columns:
        raise TypeError(
            "bool 更新要求 df 中的来源列为 bool 类型: "
            f"{non_bool_source_columns}"
        )

    # adata.obs 必须包含索引列, 否则无法和 df 逐行对齐。
    missing_obs_index_columns = [
        column for column in index_columns if column not in adata.obs
    ]
    if missing_obs_index_columns:
        raise KeyError(f"adata.obs 缺少索引列: {missing_obs_index_columns}")

    # df 是更新来源, 每个索引键只能出现一次, 否则无法确定采用哪一行的 bool 值。
    _ensure_unique_keys(df, index_columns, frame_name="df")

    # obs 重复键默认视为错误; 如果用户显式允许, 则一条 df 记录可以作用到
    # 多个 obs 行, 但保留 warning 提醒。
    duplicated_obs_keys = adata.obs.loc[:, list(index_columns)].duplicated(keep=False)
    if duplicated_obs_keys.any():
        preview = _duplicated_key_preview(adata.obs, index_columns)
        message = (
            f"adata.obs 中索引列 {list(index_columns)!r} 对应的索引不唯一; "
            f"重复示例: {preview}"
        )
        if check_adata_index_unique:
            raise ValueError(message)
        warnings.warn(message, UserWarning, stacklevel=2)

    # 根据 inplace 决定是直接修改原对象, 还是创建副本后返回。
    target = adata if inplace else adata.copy()

    # 如果目标列不存在, 先创建 NA 列, 后续只会把满足 bool 条件的行写入新值。
    missing_target_columns = [
        column for column in target_columns if column not in target.obs
    ]
    for column in missing_target_columns:
        target.obs.loc[:, column] = pd.NA
    if missing_target_columns:
        warnings.warn(
            "adata.obs 缺少目标列, 已创建并填充为 NA: "
            f"{missing_target_columns}",
            UserWarning,
            stacklevel=2,
        )

    # 用 MultiIndex 表示单列或多列组合键, 得到与 adata.obs 等长的匹配掩码。
    df_keys = pd.MultiIndex.from_frame(df.loc[:, list(index_columns)])
    obs_keys = pd.MultiIndex.from_frame(target.obs.loc[:, list(index_columns)])
    matched = obs_keys.isin(df_keys)
    if not matched.any():
        raise ValueError(
            f"adata.obs 中没有任何行能根据索引列 {list(index_columns)!r} 匹配到 df"
        )

    # 逐组处理 bool 来源列、目标列、匹配值和写入值。
    matched_positions = np.flatnonzero(matched)
    for position, (source_column, target_column) in enumerate(
        zip(source_columns, target_columns, strict=True)
    ):
        # 按 obs 匹配行顺序取出 df 中的 bool 开关; 缺失值按 False 处理。
        values_by_key = pd.Series(df[source_column].to_numpy(), index=df_keys)
        source_flags = values_by_key.reindex(obs_keys[matched]).astype(
            "boolean"
        ).fillna(False)

        # 同时要求目标列当前值等于 match_values[position]。
        target_matches = target.obs.loc[matched, target_column].eq(
            match_values[position]
        )

        # 用位置索引写入, 避免 obs index 有重复时误更新同名行。
        rows_to_update = matched_positions[
            source_flags.to_numpy(dtype=bool) & target_matches.to_numpy(dtype=bool)
        ]
        if len(rows_to_update) == 0:
            continue

        # 兼容 category 目标列: 写入前先补齐新类别。
        update_value = update_values[position]
        _ensure_assignable_categories(
            target.obs,
            target_column,
            pd.Series([update_value]),
        )
        target.obs.iloc[
            rows_to_update,
            target.obs.columns.get_loc(target_column),
        ] = update_value

    # 原位模式遵循 pandas/anndata 常见习惯返回 None; 非原位模式返回更新副本。
    if inplace:
        return None
    return target


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
