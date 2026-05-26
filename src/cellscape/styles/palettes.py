"""Reusable color palettes."""

from __future__ import annotations

from matplotlib.colors import LinearSegmentedColormap

SEURAT_HEAT_COLORS = [
    "#6c68aa",
    "#679eab",
    "#e0e9ac",
    "#fad994",
    "#f9af82",
    "#e17364",
    "#a42d51",
]
SEURAT_GREY_TO_BLUE = ["#D3D3D3", '#00008B']

def seurat_heat(
    name: str = "seurat_heat",
    *,
    n_colors: int = 256,
) -> LinearSegmentedColormap:
    """Return the continuous colormap migrated from ``raw_palettes.py``."""
    return LinearSegmentedColormap.from_list(
        name=name,
        colors=SEURAT_HEAT_COLORS,
        N=n_colors,
    )

def seurat_grey_to_blue(
    name: str = "seurat_grey_to_blue",
    *,
    n_colors: int = 256,
) -> LinearSegmentedColormap:
    """Return the continuous colormap migrated from ``raw_palettes.py``."""
    return LinearSegmentedColormap.from_list(
        name=name,
        colors=SEURAT_GREY_TO_BLUE,
        N=n_colors,
    )


def glasbey_palette(n_colors: int) -> list[str]:
    """Generate categorical colors with glasbey."""
    import glasbey

    return glasbey.create_palette(n_colors)
