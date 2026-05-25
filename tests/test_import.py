import cellscape


def test_import() -> None:
    assert isinstance(cellscape.__version__, str)
    assert callable(cellscape.spatial_scatter)
    assert callable(cellscape.umap_highlight)
