"""Basic spatial scatter example with a small in-memory AnnData object."""

import numpy as np
import pandas as pd
from anndata import AnnData

import cellscape as cs

adata = AnnData(np.array([[1.0], [2.0], [0.5]]))
adata.var_names = ["GeneA"]
adata.obs = pd.DataFrame({"cluster": ["0", "1", "0"]}, index=adata.obs_names)
adata.obsm["spatial"] = np.array([[0, 0], [1, 0], [0, 1]])

cs.spatial_scatter(adata, "cluster", show=True)
