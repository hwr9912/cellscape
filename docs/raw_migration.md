# Raw Script Migration

First-stage migration summary for the original `raw_*.py` files.

| Raw file | Original content | Category | New API |
| --- | --- | --- | --- |
| `src/cellscape/spatial/raw_expr.py` | Squidpy multi-library gene expression panel | Spatial expression plot | `cellscape.spatial_expression_panels` |
| `src/cellscape/spatial/raw_highlight_cluster.py` | Multi-panel cluster highlighting with shared legend | Spatial annotation plot | `cellscape.highlight_clusters_panels` |
| `src/cellscape/spatial/raw_hl_and_expr.py` | Gene-expression and cluster rows across library columns | Spatial expression plus annotation grid | `cellscape.highlight_and_expression_grid` |
| `src/cellscape/spatial/raw_expr_and_crop.py` | Same grid with selected-cluster crop coordinates | Spatial crop plus annotation grid | `cellscape.highlight_and_expression_grid(crop_to_selected=True)` |
| `src/cellscape/single_cell/raw_plot.py` | `umap_highlight` and Hotspot local-correlation heatmap | Single-cell visualization | `cellscape.umap_highlight`, `cellscape.local_correlation_plot` |
| `src/cellscape/styles/raw_palettes.py` | Custom continuous color map | Palette management | `cellscape.cellscape_continuous_cmap` |

The original `raw_*.py` files have been removed after migration so the package
does not ship notebook-style temporary code.
