from cellscape.spatial.plots import _try_int


def test_try_int_sorts_mixed_category_values() -> None:
    values = ["10", "2", "A", 1, "B"]

    assert sorted(values, key=_try_int) == [1, "2", "10", "A", "B"]
