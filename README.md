# cellscape

`cellscape` is a personal Python package for omics visualization tools, focused on single-cell and spatial transcriptomics workflows.

- spatial gene expression panels
- spatial cluster highlighting
- spatial expression plus annotation grids, with optional crop-to-cluster
- generic spatial scatter from `obsm["spatial"]` or `obs["x"]`/`obs["y"]`
- cell-boundary polygon plotting from `obsm["cell_border"]`
- UMAP group highlighting
- Hotspot local-correlation heatmaps
- custom continuous color palette
- labelme annotation to PNG/NPY mask conversion with label visualization

## Layout

```text
cellscape/
├── src/cellscape/
│   ├── single_cell/    # Single-cell plotting utilities
│   ├── spatial/        # Spatial transcriptomics plotting utilities
│   ├── core/           # Shared plotting primitives and validation
│   ├── io/             # Data loading and export helpers
│   ├── datasets/       # Dataset helpers and preprocessing utilities
│   ├── styles/         # Themes, palettes, and figure styling
│   └── widgets/        # Optional interactive components
├── tests/              # Unit tests
├── examples/           # Runnable examples
├── notebooks/          # Exploratory notebooks
└── docs/               # Documentation drafts
```

## Install

```bash
pip install git+https://github.com/hwr9912/cellscape.git
```

## Development

Install locally:

```bash
pip install -e .
```

Install with development tools:

```bash
pip install -e ".[dev]"
```

## Basic Usage

```python
import cellscape as cs

fig, ax = cs.spatial_scatter(adata, color="celltype", show=False)
fig, ax = cs.cell_boundary_plot(adata, color="cluster", show=False)

result = cs.labelme_to_mask(
    image_path="image.png",
    annotation_json_path="image.json",
    labels=["__ignore__", "_background_", "test"],
)

# Equivalent dataset-preprocessing entry point:
# from cellscape.datasets import labelme_to_mask, labelme_to_masks
```

All plotting functions accept `show` and `save` where applicable, and avoid
running analysis or drawing figures during package import.
