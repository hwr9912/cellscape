"""Shared plotting helpers used by multiple visualization modules."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt


def apply_publication_defaults() -> None:
    """
    指定matplotlib相关的全局图片保存参数
    """
    # # 使用等线字体
    # mpl.rcParams["font.family"] = "Arial"
    # 文字作为可编辑字符串而不是形状存储
    mpl.rcParams["pdf.fonttype"] = 42


def finish_figure(
    fig: Any,
    *,
    save: str | Path | None = None,
    show: bool = True,
    dpi: int = 300,
    bbox_inches: str = "tight",
) -> Any:
    """
    Save and optionally show a figure, then return it.
    """
    if save is not None:
        fig.savefig(save, dpi=dpi, bbox_inches=bbox_inches)
    if show:
        plt.show()
    return fig


def as_axes_array(axes: Any):
    """Return axes as a flat list-like NumPy array without importing NumPy upstream."""
    import numpy as np

    return np.asarray(axes).reshape(-1)
