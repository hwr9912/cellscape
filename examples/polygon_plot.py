"""Basic cell-boundary plotting example."""

import numpy as np
import pandas as pd
from anndata import AnnData

import cellscape as cs

adata = AnnData(np.ones((2, 1)))
adata.var_names = ["GeneA"]
adata.obs = pd.DataFrame(
    {"x": [0.0, 2.0], "y": [0.0, 1.0], "cluster": ["A", "B"]},
    index=adata.obs_names,
)
adata.obsm["cell_border"] = np.array(
    [
        [[-0.4, -0.4], [0.4, -0.4], [0.4, 0.4], [-0.4, 0.4]],
        [[-0.4, -0.4], [0.4, -0.4], [0.4, 0.4], [-0.4, 0.4]],
    ]
)

cs.cell_boundary_plot(adata, color="cluster", show=True)
