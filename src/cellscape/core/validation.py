"""Input validation helpers."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np


def require_obs_keys(adata: Any, keys: str | Sequence[str]) -> None:
    """Raise a clear error if one or more ``adata.obs`` columns are missing."""
    keys = [keys] if isinstance(keys, str) else list(keys)
    missing = [key for key in keys if key not in adata.obs]
    if missing:
        raise KeyError(f"Missing adata.obs column(s): {missing}")


def require_obsm_key(adata: Any, key: str) -> None:
    """Raise a clear error if ``adata.obsm[key]`` is missing."""
    if key not in adata.obsm:
        raise KeyError(f"Missing adata.obsm key: {key!r}")


def get_spatial_coordinates(
    adata: Any,
    *,
    spatial_key: str = "spatial",
    x_key: str | None = None,
    y_key: str | None = None,
) -> np.ndarray:
    """Return two-column spatial coordinates from ``obsm`` or ``obs`` columns."""
    if x_key is not None or y_key is not None:
        if x_key is None or y_key is None:
            raise ValueError("x_key and y_key must be provided together.")
        require_obs_keys(adata, [x_key, y_key])
        return adata.obs[[x_key, y_key]].to_numpy()

    require_obsm_key(adata, spatial_key)
    coords = np.asarray(adata.obsm[spatial_key])
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError(
            f"adata.obsm[{spatial_key!r}] must have shape (n_obs, >=2); "
            f"got {coords.shape}."
        )
    return coords[:, :2]


def require_color_key(adata: Any, key: str) -> None:
    """Validate a color key against ``obs`` columns and gene names."""
    var_names = getattr(adata, "var_names", [])
    if key not in adata.obs and key not in var_names:
        raise KeyError(
            f"{key!r} was not found in adata.obs or adata.var_names."
        )
