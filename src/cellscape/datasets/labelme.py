"""Labelme annotation preprocessing helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Iterator

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from cellscape.styles.palettes import glasbey_palette

DEFAULT_LABELS = ("__ignore__", "_background_", "test")
DEFAULT_VISUALIZATION_MAX_DIMENSION = 2048


@dataclass(frozen=True)
class LabelmeMaskResult:
    """保存单个 labelme 注释转换后的输出路径和标签映射"""

    output_dir: Path
    sample_id: str
    class_png_dir: Path
    class_png_files: dict[str, Path]
    class_npy: Path
    visualization_jpg: Path
    class_names_file: Path
    labels: tuple[str, ...]
    class_ids: dict[str, int]


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


@contextmanager
def _open_image(
    image_path: str | Path,
    *,
    allow_large_images: bool,
) -> Iterator[Image.Image]:
    """打开图片并把 Pillow 大图限制错误转换成面向用户的提示"""
    try:
        with _pillow_image_size_limit(allow_large_images):
            with Image.open(image_path) as image:
                yield image
    except Image.DecompressionBombError as exc:
        raise RuntimeError(
            "图片尺寸超过 Pillow 的安全限制"
            "如果你已经确认图片来源可信, 并且机器内存足够生成同尺寸 mask, "
            "请设置 `allow_large_images=True` 后重试"
        ) from exc


def _read_labelme_json(annotation_json_path: str | Path) -> dict[str, Any]:
    """读取并校验 labelme JSON, 至少要求存在列表类型的 shapes 字段"""
    json_path = Path(annotation_json_path)
    if not json_path.exists():
        raise FileNotFoundError(f"找不到 labelme 注释文件: {json_path}")
    with json_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if "shapes" not in data or not isinstance(data["shapes"], list):
        raise ValueError("labelme 注释文件缺少有效的 `shapes` 字段")
    return data


def _normalize_labels(labels: Sequence[str]) -> tuple[str, ...]:
    """清理标签列表, 去掉空白标签并拒绝空列表或重复标签"""
    normalized = tuple(str(label).strip() for label in labels if str(label).strip())
    if not normalized:
        raise ValueError("`labels` 至少需要包含一个有效标签")
    if len(set(normalized)) != len(normalized):
        raise ValueError("`labels` 中存在重复标签")
    return normalized


def _class_ids_from_labels(labels: Sequence[str]) -> dict[str, int]:
    """按 labelme 约定生成类别 ID, ignore=-1, background=0, 其余从 1 递增"""
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
    """从 shapes 中收集实际出现过的非空 label 名称"""
    return {str(shape["label"]) for shape in shapes if shape.get("label")}


def _draw_shape(
    draw: ImageDraw.ImageDraw,
    points: Sequence[Sequence[float]],
    shape_type: str | None,
    fill: int,
) -> None:
    """把单个 labelme shape 绘制到 Pillow ImageDraw 对象上"""
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
    """把单个 shape 栅格化成与原图同尺寸的 bool 二值 mask"""
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
    """按 labels 顺序构建每个类别的二值 mask, 并处理背景类别"""
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
    """返回单个注释文件默认的 mask 输出目录"""
    return annotation_json_path.parent / "mask"


def _write_class_names(output_dir: Path, labels: Sequence[str]) -> Path:
    """把类别名称按顺序写入 class_names.txt"""
    class_names_file = output_dir / "class_names.txt"
    class_names_file.parent.mkdir(parents=True, exist_ok=True)
    class_names_file.write_text("\n".join(labels), encoding="utf-8")
    return class_names_file


def _rgb_tuple(color: str) -> tuple[int, int, int]:
    """把 #RRGGBB 颜色字符串转换成 RGB 整数三元组"""
    color = color.lstrip("#")
    return tuple(int(color[idx : idx + 2], 16) for idx in (0, 2, 4))


def _visualization_colors(labels: Sequence[str]) -> dict[str, tuple[int, int, int]]:
    """为可视化中的每个 label 分配一个 RGB 颜色"""
    palette = glasbey_palette(len(labels))
    return {label: _rgb_tuple(color) for label, color in zip(labels, palette)}


def _scaled_size(
    size: tuple[int, int],
    *,
    max_dimension: int = DEFAULT_VISUALIZATION_MAX_DIMENSION,
) -> tuple[int, int]:
    """计算保持长宽比的可视化输出尺寸, 最长边不超过 max_dimension"""
    width, height = size
    longest = max(width, height)
    if longest <= max_dimension:
        return width, height
    scale = max_dimension / longest
    return max(1, round(width * scale)), max(1, round(height * scale))


def _load_legend_font(width: int) -> ImageFont.ImageFont:
    """按可视化图宽度选择 legend 字体大小, 并兼容旧版 Pillow"""
    font_size = max(10, min(32, round(width * 0.015)))
    try:
        return ImageFont.load_default(size=font_size)
    except TypeError:
        return ImageFont.load_default()


def _resized_mask(mask: np.ndarray, size: tuple[int, int]) -> np.ndarray:
    """把原尺寸 bool mask 用最近邻缩放到可视化尺寸"""
    mask_img = Image.fromarray(mask.astype(np.uint8) * 255)
    if mask_img.size != size:
        mask_img = mask_img.resize(size, Image.Resampling.NEAREST)
    return np.asarray(mask_img, dtype=np.uint8) > 0


def _legend_metrics(
    draw: ImageDraw.ImageDraw,
    labels: Sequence[str],
    font: ImageFont.ImageFont,
    image_width: int,
) -> tuple[int, int, int, int, int, int]:
    """计算 legend 盒子和内部间距, 宽度约束在图宽 5-10%"""
    target_width = max(1, round(image_width * 0.075))
    padding = max(4, round(target_width * 0.10))
    swatch = max(6, round(target_width * 0.16))
    text_gap = max(3, round(target_width * 0.08))
    row_gap = max(3, round(target_width * 0.06))
    text_sizes = [draw.textbbox((0, 0), label, font=font) for label in labels]
    text_width = max((bbox[2] - bbox[0] for bbox in text_sizes), default=0)
    text_height = max((bbox[3] - bbox[1] for bbox in text_sizes), default=10)
    legend_width = padding * 2 + swatch + text_gap + text_width
    min_width = max(1, round(image_width * 0.05))
    max_width = max(min_width, round(image_width * 0.10))
    legend_width = min(max(legend_width, min_width), max_width)
    return legend_width, padding, swatch, text_gap, row_gap, text_height


def _save_visualization(
    image_path: Path,
    masks: Mapping[str, np.ndarray],
    labels: Sequence[str],
    output_path: Path,
    *,
    alpha: float,
    allow_large_images: bool,
) -> None:
    """保存带 mask 叠加和 legend 的预览 JPG, 大图会先等比例缩小"""
    with _pillow_image_size_limit(allow_large_images):
        with _open_image(
            image_path,
            allow_large_images=allow_large_images,
        ) as image:
            target_size = _scaled_size(image.size)
            if image.size != target_size:
                image.thumbnail(target_size, Image.Resampling.LANCZOS)
            base = image.convert("RGBA")
        overlay = Image.new("RGBA", base.size, (0, 0, 0, 0))
        overlay_array = np.asarray(overlay).copy()
        colors = _visualization_colors(labels)

        for label in labels:
            r, g, b = colors[label]
            overlay_array[_resized_mask(masks[label], base.size)] = (
                r,
                g,
                b,
                int(round(alpha * 255)),
            )

        blended = Image.alpha_composite(base, Image.fromarray(overlay_array, mode="RGBA"))
        draw = ImageDraw.Draw(blended)
        font = _load_legend_font(blended.width)
        (
            legend_width,
            padding,
            swatch,
            text_gap,
            row_gap,
            text_height,
        ) = _legend_metrics(draw, labels, font, blended.width)
        legend_height = (
            padding * 2
            + len(labels) * text_height
            + max(len(labels) - 1, 0) * row_gap
        )
        margin = max(4, round(blended.width * 0.008))
        x0 = max(0, blended.width - legend_width - margin)
        y0 = max(0, blended.height - legend_height - margin)
        x1 = max(x0, blended.width - margin)
        y1 = max(y0, blended.height - margin)

        draw.rounded_rectangle(
            [x0, y0, x1, y1],
            radius=max(2, round(blended.width * 0.002)),
            fill=(255, 255, 255, 210),
        )
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
    """解析原图路径, 优先使用显式 image_path, 否则使用 JSON 中的 imagePath"""
    if image_path is not None:
        resolved = Path(image_path)
    else:
        image_name = data.get("imagePath")
        if not image_name:
            raise ValueError(f"{json_path} 缺少 `imagePath`, 请显式传入 `image_path`")
        resolved = json_path.parent / str(image_name)
    if not resolved.exists():
        raise FileNotFoundError(f"找不到输入图片: {resolved}")
    return resolved


def _sample_id_from_json(annotation_json_path: Path, sample_id: str | None) -> str:
    """返回输出样本 ID; 未显式传入时使用 JSON 文件名 stem"""
    return sample_id or annotation_json_path.stem


def labelme_to_mask(
    image_path: str | Path | None,
    annotation_json_path: str | Path,
    *,
    labels: Sequence[str] = DEFAULT_LABELS,
    output_dir: str | Path | None = None,
    sample_id: str | None = None,
    visualization_alpha: float | None = 0.35,
    allow_large_images: bool = False,
) -> LabelmeMaskResult:
    """
    将单个 labelme 注释文件转换为按类别拆分的灰度二值 mask 数据集

    PNG 使用 0/255 的 uint8 编码, 每个 label 一个文件; NPY 使用 0/1 的 uint8 编码, 
    每个样本一个文件, 按 `labels` 顺序在不同维度存储

    ### Parameters
    - `image_path: str | Path | None` 原始图片路径; 为 None 时使用 json 内的 `imagePath`,
    - `annotation_json_path: str | Path` labelme 生成的 json 注释文件路径,
    - `labels: Sequence[str] = DEFAULT_LABELS` label 顺序, 遵循官方 `__ignore__`、`_background_` 约定,
    - `output_dir: str | Path | None = None` 输出根目录; 默认是 json 同目录下的 `mask`,
    - `sample_id: str | None = None` 输出样本名; 默认使用 json 文件名,
    - `visualization_alpha: float | None = 0.35` 可视化叠加透明度; 为 None 时跳过可视化,
    - `allow_large_images: bool = False` 是否忽略 Pillow 图片大小安全限制; 
      仅在确认图片可信且内存足够生成同尺寸 mask 时设置为 True

    ### Example
    ```python
    import cellscape.datasets as dts

    dts.labelme_to_mask(
        image_path="data/labelme_test/random1.tif",
        annotation_json_path="data/labelme_test/random1.json",
        labels=["__ignore__", "_background_", "test"],
        output_dir="data/labelme_test/mask",
    )
    ```
    """
    annotation_json_path = Path(annotation_json_path)
    output_dir = Path(output_dir) if output_dir is not None else _default_output_dir(annotation_json_path)
    labels = _normalize_labels(labels)
    class_ids = _class_ids_from_labels(labels)

    data = _read_labelme_json(annotation_json_path)
    image_path = _resolve_image_path(annotation_json_path, data, image_path)
    with _open_image(image_path, allow_large_images=allow_large_images) as img:
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
    if visualization_alpha is not None:
        _save_visualization(
            image_path,
            masks,
            labels,
            visualization_jpg,
            alpha=visualization_alpha,
            allow_large_images=allow_large_images,
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
    visualization_alpha: float | None = 0.35,
    allow_large_images: bool = False,
) -> list[LabelmeMaskResult]:
    """
    批量将 labelme 注释 json 转换为按类别拆分的灰度二值 mask 数据集

    每个 json 默认使用其中的 `imagePath` 查找图片; 如果传入 `image_dir`, 则优先在
    该目录下按 `imagePath` 的文件名匹配输出根目录默认是注释目录下的 `mask`

    ### Parameters
    - `annotation_dir: str | Path` labelme json 文件所在目录,
    - `image_dir: str | Path | None = None` 图片目录; 为 None 时使用 json 内的 `imagePath`,
    - `output_dir: str | Path | None = None` 输出根目录; 为 None 时输出到 annotation_dir/mask,
    - `labels: Sequence[str] = DEFAULT_LABELS` label 顺序, 遵循官方 `__ignore__`、`_background_` 约定,
    - `pattern: str = "*.json"` json 文件匹配模式,
    - `recursive: bool = False` 是否递归扫描 annotation_dir,
    - `visualization_alpha: float | None = 0.35` 可视化叠加透明度; 为 None 时跳过可视化,
    - `allow_large_images: bool = False` 是否忽略 Pillow 图片大小安全限制; 
      仅在确认图片可信且内存足够生成同尺寸 mask 时设置为 True

    ### Example
    ```python
    import cellscape.datasets as dts

    dts.labelme_to_masks(
        annotation_dir="library_backgrounds",
        image_dir="library_backgrounds",
        output_dir="library_backgrounds/mask",
        labels=["__ignore__", "_background_", "test"],
    )
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
                raise ValueError(f"{json_path} 缺少 `imagePath`, 无法在 image_dir 中匹配图片")
            resolved_image_path = Path(image_dir) / Path(str(image_name)).name

        results.append(
            labelme_to_mask(
                resolved_image_path,
                json_path,
                labels=labels,
                output_dir=output_dir,
                visualization_alpha=visualization_alpha,
                allow_large_images=allow_large_images,
            )
        )

    return results


batch_labelme_to_masks = labelme_to_masks
