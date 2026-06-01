"""Labelme annotation preprocessing helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from cellscape.styles.palettes import glasbey_palette

DEFAULT_LABELS = ("__ignore__", "_background_", "test")


@dataclass(frozen=True)
class LabelmeMaskResult:
    """保存单个 labelme 注释转换后的输出路径和标签映射。"""

    output_dir: Path
    sample_id: str
    class_png_dir: Path
    class_png_files: dict[str, Path]
    class_npy: Path
    visualization_jpg: Path
    class_names_file: Path
    labels: tuple[str, ...]
    class_ids: dict[str, int]


def _read_labelme_json(annotation_json_path: str | Path) -> dict[str, Any]:
    json_path = Path(annotation_json_path)
    if not json_path.exists():
        raise FileNotFoundError(f"找不到 labelme 注释文件: {json_path}")
    with json_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if "shapes" not in data or not isinstance(data["shapes"], list):
        raise ValueError("labelme 注释文件缺少有效的 `shapes` 字段。")
    return data


def _normalize_labels(labels: Sequence[str]) -> tuple[str, ...]:
    normalized = tuple(str(label).strip() for label in labels if str(label).strip())
    if not normalized:
        raise ValueError("`labels` 至少需要包含一个有效标签。")
    if len(set(normalized)) != len(normalized):
        raise ValueError("`labels` 中存在重复标签。")
    return normalized


def _class_ids_from_labels(labels: Sequence[str]) -> dict[str, int]:
    normalized = _normalize_labels(labels)
    class_ids: dict[str, int] = {}
    next_positive_id = 1
    for label in normalized:
        if label == "__ignore__":
            class_ids[label] = -1
        elif label == "_background_":
            class_ids[label] = 0
        else:
            class_ids[label] = next_positive_id
            next_positive_id += 1
    return class_ids


def _labels_from_shapes(shapes: Sequence[Mapping[str, Any]]) -> set[str]:
    return {str(shape["label"]) for shape in shapes if shape.get("label")}


def _draw_shape(
    draw: ImageDraw.ImageDraw,
    points: Sequence[Sequence[float]],
    shape_type: str | None,
    fill: int,
) -> None:
    xy = [tuple(map(float, point[:2])) for point in points]
    if not xy:
        return

    if shape_type == "rectangle" and len(xy) >= 2:
        x_values = [point[0] for point in xy[:2]]
        y_values = [point[1] for point in xy[:2]]
        draw.rectangle(
            [min(x_values), min(y_values), max(x_values), max(y_values)],
            fill=fill,
        )
    elif shape_type == "circle" and len(xy) >= 2:
        center, edge = xy[0], xy[1]
        radius = float(np.hypot(edge[0] - center[0], edge[1] - center[1]))
        draw.ellipse(
            [
                center[0] - radius,
                center[1] - radius,
                center[0] + radius,
                center[1] + radius,
            ],
            fill=fill,
        )
    elif shape_type in {"line", "linestrip"} and len(xy) >= 2:
        draw.line(xy, fill=fill, width=1)
    elif shape_type == "point":
        x, y = xy[0]
        draw.ellipse([x - 1, y - 1, x + 1, y + 1], fill=fill)
    elif len(xy) >= 3:
        draw.polygon(xy, fill=fill)


def _shape_to_binary_mask(
    width: int,
    height: int,
    shape: Mapping[str, Any],
) -> np.ndarray:
    points = shape.get("points")
    if not isinstance(points, list):
        return np.zeros((height, width), dtype=bool)
    mask_img = Image.new("L", (width, height), 0)
    draw = ImageDraw.Draw(mask_img)
    _draw_shape(draw, points, shape.get("shape_type"), 1)
    return np.asarray(mask_img, dtype=bool)


def _build_binary_masks(
    shapes: Sequence[Mapping[str, Any]],
    labels: Sequence[str],
    width: int,
    height: int,
) -> dict[str, np.ndarray]:
    shape_labels = _labels_from_shapes(shapes)
    unknown_labels = sorted(shape_labels.difference(labels))
    if unknown_labels:
        raise ValueError(f"json 中存在未写入 `labels` 的标签: {unknown_labels}")

    masks = {
        label: np.zeros((height, width), dtype=bool)
        for label in labels
    }

    for shape in shapes:
        label = shape.get("label")
        if not label:
            continue
        label_str = str(label)
        if label_str not in masks:
            continue
        masks[label_str] |= _shape_to_binary_mask(width, height, shape)

    positive_union = np.zeros((height, width), dtype=bool)
    for label, mask in masks.items():
        if label not in {"__ignore__", "_background_"}:
            positive_union |= mask

    if "__ignore__" in masks and "__ignore__" not in shape_labels:
        masks["__ignore__"] = np.zeros((height, width), dtype=bool)
    if "_background_" in masks and "_background_" not in shape_labels:
        masks["_background_"] = ~positive_union

    return masks


def _default_output_dir(annotation_json_path: Path) -> Path:
    return annotation_json_path.parent / "mask"


def _write_class_names(output_dir: Path, labels: Sequence[str]) -> Path:
    class_names_file = output_dir / "class_names.txt"
    class_names_file.parent.mkdir(parents=True, exist_ok=True)
    class_names_file.write_text("\n".join(labels), encoding="utf-8")
    return class_names_file


def _rgb_tuple(color: str) -> tuple[int, int, int]:
    color = color.lstrip("#")
    return tuple(int(color[idx : idx + 2], 16) for idx in (0, 2, 4))


def _visualization_colors(labels: Sequence[str]) -> dict[str, tuple[int, int, int]]:
    palette = glasbey_palette(len(labels))
    return {label: _rgb_tuple(color) for label, color in zip(labels, palette)}


def _save_visualization(
    image_path: Path,
    masks: Mapping[str, np.ndarray],
    labels: Sequence[str],
    output_path: Path,
    *,
    alpha: float,
) -> None:
    base = Image.open(image_path).convert("RGBA")
    overlay = Image.new("RGBA", base.size, (0, 0, 0, 0))
    overlay_array = np.asarray(overlay).copy()
    colors = _visualization_colors(labels)

    for label in labels:
        r, g, b = colors[label]
        overlay_array[masks[label]] = (r, g, b, int(round(alpha * 255)))

    blended = Image.alpha_composite(base, Image.fromarray(overlay_array, mode="RGBA"))
    draw = ImageDraw.Draw(blended)
    font = ImageFont.load_default()

    swatch = 12
    padding = 8
    row_gap = 5
    text_gap = 6
    text_sizes = [draw.textbbox((0, 0), label, font=font) for label in labels]
    text_width = max((bbox[2] - bbox[0] for bbox in text_sizes), default=0)
    text_height = max((bbox[3] - bbox[1] for bbox in text_sizes), default=10)
    legend_width = padding * 2 + swatch + text_gap + text_width
    legend_height = padding * 2 + len(labels) * text_height + max(len(labels) - 1, 0) * row_gap
    x0 = max(0, blended.width - legend_width - 8)
    y0 = max(0, blended.height - legend_height - 8)
    x1 = blended.width - 8
    y1 = blended.height - 8

    draw.rounded_rectangle([x0, y0, x1, y1], radius=4, fill=(255, 255, 255, 210))
    for row_idx, label in enumerate(labels):
        y = y0 + padding + row_idx * (text_height + row_gap)
        draw.rectangle(
            [x0 + padding, y, x0 + padding + swatch, y + swatch],
            fill=colors[label],
        )
        draw.text(
            (x0 + padding + swatch + text_gap, y - 1),
            label,
            fill=(0, 0, 0, 255),
            font=font,
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    blended.convert("RGB").save(output_path, quality=95)


def _resolve_image_path(
    json_path: Path,
    data: Mapping[str, Any],
    image_path: str | Path | None,
) -> Path:
    if image_path is not None:
        resolved = Path(image_path)
    else:
        image_name = data.get("imagePath")
        if not image_name:
            raise ValueError(f"{json_path} 缺少 `imagePath`，请显式传入 `image_path`。")
        resolved = json_path.parent / str(image_name)
    if not resolved.exists():
        raise FileNotFoundError(f"找不到输入图片: {resolved}")
    return resolved


def _sample_id_from_json(annotation_json_path: Path, sample_id: str | None) -> str:
    return sample_id or annotation_json_path.stem


def labelme_to_mask(
    image_path: str | Path | None,
    annotation_json_path: str | Path,
    *,
    labels: Sequence[str] = DEFAULT_LABELS,
    output_dir: str | Path | None = None,
    sample_id: str | None = None,
    visualization_alpha: float = 0.35,
) -> LabelmeMaskResult:
    """
    将单个 labelme 注释文件转换为按类别拆分的灰度二值 mask 数据集

    PNG 使用 0/255 的 uint8 编码，每个 label 一个文件；NPY 使用 0/1 的 uint8 编码，
    每个样本一个文件，按 `labels` 顺序在不同维度存储。

    ### Parameters
    - `image_path: str | Path | None` 原始图片路径；为 None 时使用 json 内的 `imagePath`,
    - `annotation_json_path: str | Path` labelme 生成的 json 注释文件路径,
    - `labels: Sequence[str] = DEFAULT_LABELS` label 顺序，遵循官方 `__ignore__`、`_background_` 约定,
    - `output_dir: str | Path | None = None` 输出根目录；默认是 json 同目录下的 `mask`,
    - `sample_id: str | None = None` 输出样本名；默认使用 json 文件名,
    - `visualization_alpha: float = 0.35` 可视化叠加透明度

    ### Example
    ```python

    ```
    """
    annotation_json_path = Path(annotation_json_path)
    output_dir = Path(output_dir) if output_dir is not None else _default_output_dir(annotation_json_path)
    labels = _normalize_labels(labels)
    class_ids = _class_ids_from_labels(labels)

    data = _read_labelme_json(annotation_json_path)
    image_path = _resolve_image_path(annotation_json_path, data, image_path)
    with Image.open(image_path) as img:
        width, height = img.size

    sample_id = _sample_id_from_json(annotation_json_path, sample_id)
    masks = _build_binary_masks(data["shapes"], labels, width, height)

    class_png_dir = output_dir / "SegmentationClass" / sample_id
    class_npy = output_dir / "SegmentationClassNpy" / f"{sample_id}.npy"
    visualization_jpg = output_dir / "SegmentationClassVisualization" / f"{sample_id}.jpg"
    class_names_file = _write_class_names(output_dir, labels)

    class_png_dir.mkdir(parents=True, exist_ok=True)
    class_png_files: dict[str, Path] = {}
    stacked_masks = []
    for label in labels:
        class_id = class_ids[label]
        png_path = class_png_dir / f"{class_id}_{label}.png"
        Image.fromarray((masks[label].astype(np.uint8) * 255)).save(png_path)
        class_png_files[label] = png_path
        stacked_masks.append(masks[label].astype(np.uint8))

    class_npy.parent.mkdir(parents=True, exist_ok=True)
    np.save(class_npy, np.stack(stacked_masks, axis=0).astype(np.uint8))
    _save_visualization(
        image_path,
        masks,
        labels,
        visualization_jpg,
        alpha=visualization_alpha,
    )

    return LabelmeMaskResult(
        output_dir=output_dir,
        sample_id=sample_id,
        class_png_dir=class_png_dir,
        class_png_files=class_png_files,
        class_npy=class_npy,
        visualization_jpg=visualization_jpg,
        class_names_file=class_names_file,
        labels=labels,
        class_ids=class_ids,
    )


def labelme_to_masks(
    annotation_dir: str | Path,
    *,
    image_dir: str | Path | None = None,
    output_dir: str | Path | None = None,
    labels: Sequence[str] = DEFAULT_LABELS,
    pattern: str = "*.json",
    recursive: bool = False,
    visualization_alpha: float = 0.35,
) -> list[LabelmeMaskResult]:
    """
    批量将 labelme 注释 json 转换为按类别拆分的灰度二值 mask 数据集

    每个 json 默认使用其中的 `imagePath` 查找图片；如果传入 `image_dir`，则优先在
    该目录下按 `imagePath` 的文件名匹配。输出根目录默认是注释目录下的 `mask`。

    ### Parameters
    - `annotation_dir: str | Path` labelme json 文件所在目录,
    - `image_dir: str | Path | None = None` 图片目录；为 None 时使用 json 内的 `imagePath`,
    - `output_dir: str | Path | None = None` 输出根目录；为 None 时输出到 annotation_dir/mask,
    - `labels: Sequence[str] = DEFAULT_LABELS` label 顺序，遵循官方 `__ignore__`、`_background_` 约定,
    - `pattern: str = "*.json"` json 文件匹配模式,
    - `recursive: bool = False` 是否递归扫描 annotation_dir,
    - `visualization_alpha: float = 0.35` 可视化叠加透明度

    ### Example
    ```python

    ```
    """
    annotation_dir = Path(annotation_dir)
    if not annotation_dir.exists():
        raise FileNotFoundError(f"找不到注释目录: {annotation_dir}")

    json_paths = sorted(annotation_dir.rglob(pattern) if recursive else annotation_dir.glob(pattern))
    if not json_paths:
        raise ValueError(f"在目录中没有找到匹配 `{pattern}` 的 json 文件: {annotation_dir}")

    output_dir = Path(output_dir) if output_dir is not None else annotation_dir / "mask"
    labels = _normalize_labels(labels)
    _write_class_names(output_dir, labels)

    results: list[LabelmeMaskResult] = []
    for json_path in json_paths:
        data = _read_labelme_json(json_path)
        resolved_image_path = None
        if image_dir is not None:
            image_name = data.get("imagePath")
            if not image_name:
                raise ValueError(f"{json_path} 缺少 `imagePath`，无法在 image_dir 中匹配图片。")
            resolved_image_path = Path(image_dir) / Path(str(image_name)).name

        results.append(
            labelme_to_mask(
                resolved_image_path,
                json_path,
                labels=labels,
                output_dir=output_dir,
                visualization_alpha=visualization_alpha,
            )
        )

    return results


batch_labelme_to_masks = labelme_to_masks
