from cellscape.spatial.plots import _try_int
from cellscape.spatial.annotation import _load_png_masks
from PIL import Image


def test_try_int_sorts_mixed_category_values() -> None:
    values = ["10", "2", "A", 1, "B"]

    assert sorted(values, key=_try_int) == [1, "2", "10", "A", "B"]


def test_load_png_masks_allows_large_trusted_masks(tmp_path) -> None:
    previous_limit = Image.MAX_IMAGE_PIXELS
    Image.MAX_IMAGE_PIXELS = 40
    try:
        Image.new("L", (10, 10), 255).save(tmp_path / "region.png")

        masks = _load_png_masks(tmp_path, ["region"])
    finally:
        Image.MAX_IMAGE_PIXELS = previous_limit

    assert masks.shape == (1, 10, 10)
    assert masks.dtype == bool
    assert masks.all()
