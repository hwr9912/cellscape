"""Basic UMAP highlight example.

This example assumes ``adata`` already contains ``obsm['X_umap']`` and the
requested obs columns.
"""

import matplotlib.pyplot as plt

import cellscape as cs


def plot_group(adata):
    fig, ax = plt.subplots(figsize=(4, 4))
    cs.umap_highlight(
        adata,
        ax,
        group_value="treated",
        group_col="condition",
        category_col="celltype",
    )
    return fig, ax
