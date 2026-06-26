from cellscape.spatial.plots import _try_int
from cellscape.spatial.plots import highlight_clusters_panels
from cellscape.spatial.annotation import _load_png_masks
from cellscape.spatial.annotation import update_obs_from_df
from cellscape.spatial.annotation import update_obs_from_bool_df
from PIL import Image
import pandas as pd
import pytest


class _MiniAnnData:
    obs = pd.DataFrame({"celltype": ["A", "B"]})
    var_names = ["Ccl3"]

    def __init__(self, obs: pd.DataFrame | None = None) -> None:
        if obs is not None:
            self.obs = obs

    def copy(self) -> "_MiniAnnData":
        return _MiniAnnData(self.obs.copy())


def test_try_int_sorts_mixed_category_values() -> None:
    """混合数字和字符串类别时, 数字值应优先按数值顺序排序。"""
    values = ["10", "2", "A", 1, "B"]

    assert sorted(values, key=_try_int) == [1, "2", "10", "A", "B"]


def test_highlight_clusters_panels_rejects_gene_color_key() -> None:
    """聚类高亮面板不应把基因名当作 obs 色彩分组列。"""
    with pytest.raises(KeyError, match="基因表达量面板"):
        highlight_clusters_panels(
            _MiniAnnData(),
            color_key="Ccl3",
            select_cluster=None,
            panels=["library0"],
            show=False,
        )


def test_load_png_masks_allows_large_trusted_masks(tmp_path) -> None:
    """读取可信 mask PNG 时应临时放宽 Pillow 的像素安全限制。"""
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


def test_update_obs_from_df_updates_matched_rows_only() -> None:
    """只更新能按索引列匹配到 df 的 obs 行, 未匹配行保持原值。"""
    adata = _MiniAnnData(
        pd.DataFrame(
            {
                "library": ["lib0", "lib0", "lib1"],
                "cell_id": ["c1", "c2", "c3"],
                "x": [1.0, 2.0, 3.0],
                "y": [4.0, 5.0, 6.0],
            }
        )
    )
    df = pd.DataFrame(
        {
            "library": ["lib0", "lib1"],
            "cell_id": ["c2", "c3"],
            "new_x": [20.0, 30.0],
            "new_y": [50.0, 60.0],
        }
    )

    update_obs_from_df(
        adata,
        df,
        index_columns=["library", "cell_id"],
        source_columns=["new_x", "new_y"],
        target_columns=["x", "y"],
    )

    assert adata.obs["x"].tolist() == [1.0, 20.0, 30.0]
    assert adata.obs["y"].tolist() == [4.0, 50.0, 60.0]


def test_update_obs_from_df_creates_missing_target_column() -> None:
    """目标列不存在时应创建 NA 列并仅写入匹配行的新值。"""
    adata = _MiniAnnData(pd.DataFrame({"cell_id": ["c1", "c2"]}))
    df = pd.DataFrame({"cell_id": ["c2"], "new_region": ["tumor"]})

    with pytest.warns(UserWarning, match="缺少目标列"):
        update_obs_from_df(
            adata,
            df,
            index_columns="cell_id",
            source_columns="new_region",
            target_columns="region",
        )

    assert pd.isna(adata.obs.loc[0, "region"])
    assert adata.obs.loc[1, "region"] == "tumor"


def test_update_obs_from_df_rejects_duplicate_df_keys() -> None:
    """df 的索引列组合必须唯一, 否则无法确定用于更新的来源行。"""
    adata = _MiniAnnData(pd.DataFrame({"cell_id": ["c1"], "region": ["old"]}))
    df = pd.DataFrame({"cell_id": ["c1", "c1"], "region": ["A", "B"]})

    with pytest.raises(ValueError, match="df.*索引不唯一"):
        update_obs_from_df(
            adata,
            df,
            index_columns="cell_id",
            source_columns="region",
            target_columns="region",
        )


def test_update_obs_from_df_rejects_empty_match() -> None:
    """没有任何 obs 行能匹配到 df 时应报错, 避免静默跳过更新。"""
    adata = _MiniAnnData(pd.DataFrame({"cell_id": ["c1"], "region": ["old"]}))
    df = pd.DataFrame({"cell_id": ["c2"], "region": ["new"]})

    with pytest.raises(ValueError, match="没有任何行.*匹配到 df"):
        update_obs_from_df(
            adata,
            df,
            index_columns="cell_id",
            source_columns="region",
            target_columns="region",
        )


def test_update_obs_from_df_warns_duplicate_obs_keys_when_allowed() -> None:
    """允许 obs 重复索引时仍应发出 warning, 并更新所有匹配行。"""
    adata = _MiniAnnData(
        pd.DataFrame({"cell_id": ["c1", "c1"], "region": ["old", "old"]})
    )
    df = pd.DataFrame({"cell_id": ["c1"], "region": ["new"]})

    with pytest.warns(UserWarning, match="adata.obs.*索引不唯一"):
        update_obs_from_df(
            adata,
            df,
            index_columns="cell_id",
            source_columns="region",
            target_columns="region",
            check_adata_index_unique=False,
        )

    assert adata.obs["region"].tolist() == ["new", "new"]


def test_update_obs_from_df_can_return_copy() -> None:
    """非原位模式应返回更新后的副本, 并保持原始 adata 不变。"""
    adata = _MiniAnnData(pd.DataFrame({"cell_id": ["c1"], "region": ["old"]}))
    df = pd.DataFrame({"cell_id": ["c1"], "region": ["new"]})

    result = update_obs_from_df(
        adata,
        df,
        index_columns="cell_id",
        source_columns="region",
        target_columns="region",
        inplace=False,
    )

    assert adata.obs.loc[0, "region"] == "old"
    assert result.obs.loc[0, "region"] == "new"


def test_update_obs_from_bool_df_updates_true_matches_only() -> None:
    """bool 更新只在来源列为 True 且目标列当前值匹配时写入新值。"""
    adata = _MiniAnnData(
        pd.DataFrame(
            {
                "cell_id": ["c1", "c2", "c3", "c4"],
                "region": ["old", "old", "keep", "old"],
                "state": ["low", "low", "low", "high"],
            }
        )
    )
    df = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "region_flag": [True, False, True, True],
            "state_flag": [False, True, True, True],
        }
    )

    update_obs_from_bool_df(
        adata,
        df,
        index_columns="cell_id",
        source_columns=["region_flag", "state_flag"],
        target_columns=["region", "state"],
        match_values=["old", "low"],
        update_values=["new", "high"],
    )

    assert adata.obs["region"].tolist() == ["new", "old", "keep", "new"]
    assert adata.obs["state"].tolist() == ["low", "high", "high", "high"]


def test_update_obs_from_bool_df_requires_bool_source() -> None:
    """bool 更新要求 df 来源列本身是 bool 类型。"""
    adata = _MiniAnnData(pd.DataFrame({"cell_id": ["c1"], "region": ["old"]}))
    df = pd.DataFrame({"cell_id": ["c1"], "region_flag": ["yes"]})

    with pytest.raises(TypeError, match="bool.*来源列为 bool 类型"):
        update_obs_from_bool_df(
            adata,
            df,
            index_columns="cell_id",
            source_columns="region_flag",
            target_columns="region",
            match_values="old",
            update_values="new",
        )


def test_update_obs_from_bool_df_requires_matching_shapes() -> None:
    """来源列、目标列以及 bool 参数必须同时为字符串或等长列表。"""
    adata = _MiniAnnData(pd.DataFrame({"cell_id": ["c1"], "region": ["old"]}))
    df = pd.DataFrame({"cell_id": ["c1"], "region_flag": [True]})

    with pytest.raises(ValueError, match="source_columns.*target_columns"):
        update_obs_from_bool_df(
            adata,
            df,
            index_columns="cell_id",
            source_columns="region_flag",
            target_columns=["region"],
            match_values="old",
            update_values="new",
        )

    with pytest.raises(ValueError, match="bool 更新.*等长列表"):
        update_obs_from_bool_df(
            adata,
            df,
            index_columns="cell_id",
            source_columns=["region_flag"],
            target_columns=["region"],
            match_values="old",
            update_values="new",
        )
