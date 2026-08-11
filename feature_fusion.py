from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd


GSE_IDS = [
    "GSE25066", "GSE194040", "GSE16716", "GSE41998", "GSE180962",
    "GSE20271", "GSE34138", "GSE50948", "GSE149322", "GSE22226",
    "GSE32603", "GSE22358", "GSE32646", "GSE16446", "GSE130788",
    "GSE231629", "GSE123845", "GSE4779", "GSE173839", "GSE22093",
    "GSE21997", "GSE42822", "GSE66399", "GSE181574", "GSE23988",
    "GSE41656", "GSE8465", "GSE122630", "GSE21974", "GSE207248",
    "GSE191127",
]

CLINICAL_SCHEMA = {
    "age":         {"type": "numeric"},
    "stage":       {"type": "categorical", "categories": ["I", "II", "III"]},
    "er_status":   {"type": "categorical", "categories": ["positive", "negative"]},
    "pr_status":   {"type": "categorical", "categories": ["positive", "negative"]},
    "her2_status": {"type": "categorical", "categories": ["positive", "negative"]},
    "grade":       {"type": "categorical", "categories": ["1", "2", "3"]},
    "gse_id":      {"type": "categorical", "categories": GSE_IDS},
}

_NA_STRINGS = {"NA", "nan", "None", ""}


def _clean_categorical(series: pd.Series, strip_trailing_zero: bool = False) -> pd.Series:
    """str.strip() only -- deliberately NOT .lower(), so values that are
    genuinely case-distinct in the schema's category list still match.
    Also normalises the string spellings of missing data to np.nan.
    """
    s = series.astype(str).str.strip()
    if strip_trailing_zero:
        
        s = s.str.replace(r"\.0$", "", regex=True)
    s = s.where(~s.isin(_NA_STRINGS), other=np.nan)
    return s


def _encode_clinical(clinical_df: pd.DataFrame) -> pd.DataFrame:
    """One-hot encode clinical_df per CLINICAL_SCHEMA. Every categorical
    column always produces one dummy column per FIXED category, regardless
    of which categories are actually observed in clinical_df -- this is
    what keeps train/test column sets identical."""
    parts = []

    for col, spec in CLINICAL_SCHEMA.items():
        if col not in clinical_df.columns:
            warnings.warn(f"[CLINICAL] Column '{col}' not found; treating as all-missing.")
            source = pd.Series(np.nan, index=clinical_df.index)
        else:
            source = clinical_df[col]

        if spec["type"] == "numeric":
            parts.append(pd.to_numeric(source, errors="coerce").rename(col))
            continue

        categories = spec["categories"]
        cleaned = _clean_categorical(source, strip_trailing_zero=(col == "grade"))
        cat_series = pd.Categorical(cleaned, categories=categories)
        dummies = pd.get_dummies(cat_series, prefix=col)

        expected_cols = [f"{col}_{c}" for c in categories]
        for c in expected_cols:
            if c not in dummies.columns:
                dummies[c] = 0
        dummies = dummies[expected_cols]
        dummies.index = clinical_df.index

        parts.append(dummies)

    return pd.concat(parts, axis=1)


def load_hybrid_scores(scores_path: str, le_path: str, tf_path: str) -> pd.DataFrame:
    """Load the pathway ssGSEA, leading-edge gene, and TF activity score
    CSVs for one cohort and join them on sample_id."""
    scores_df = pd.read_csv(scores_path, index_col=0)
    le_df     = pd.read_csv(le_path, index_col=0)
    tf_df     = pd.read_csv(tf_path, index_col=0)

    for df in (scores_df, le_df, tf_df):
        df.index = df.index.astype(str)
        df.index.name = "sample_id"

    id_sets = [set(scores_df.index), set(le_df.index), set(tf_df.index)]
    common  = set.intersection(*id_sets)
    union   = set.union(*id_sets)

    if len(common) < len(union):
        missing = sorted(union - common)
        warnings.warn(
            f"[HYBRID] {len(missing)} sample(s) not present in all three score files "
            f"(pathway/leading-edge/TF); dropping to the common set. First few: {missing[:10]}"
        )

    common_sorted = sorted(common)
    hybrid = (
        scores_df.loc[common_sorted]
        .join(le_df.loc[common_sorted], how="inner", rsuffix="_le")
        .join(tf_df.loc[common_sorted], how="inner", rsuffix="_tf")
    )
    hybrid.index.name = "sample_id"
    return hybrid


def fuse_features(clinical_path: str, scores_path: str, le_path: str,
                   tf_path: str, cohort_label: str) -> pd.DataFrame:
    """Load + one-hot-encode clinical data, load the transcriptomic score
    blocks, and inner-join them on sample_id. Only patients present on both
    sides survive; drops on either side are warned about."""

    clinical_raw = pd.read_csv(clinical_path)
    if "sample_id" not in clinical_raw.columns:
        raise ValueError(f"[FUSION] '{clinical_path}' has no 'sample_id' column.")
    clinical_raw = clinical_raw.set_index("sample_id")
    clinical_raw.index = clinical_raw.index.astype(str)

    clinical_encoded = _encode_clinical(clinical_raw)
    hybrid_scores     = load_hybrid_scores(scores_path, le_path, tf_path)

    common_ids           = clinical_encoded.index.intersection(hybrid_scores.index)
    clinical_only_ids    = clinical_encoded.index.difference(hybrid_scores.index)
    transcriptomic_only  = hybrid_scores.index.difference(clinical_encoded.index)

    if len(clinical_only_ids) > 0:
        warnings.warn(
            f"[FUSION][{cohort_label}] {len(clinical_only_ids)} sample(s) had clinical data "
            f"but no transcriptomic scores; dropped. First few: {sorted(clinical_only_ids)[:10]}"
        )
    if len(transcriptomic_only) > 0:
        warnings.warn(
            f"[FUSION][{cohort_label}] {len(transcriptomic_only)} sample(s) had transcriptomic "
            f"scores but no clinical data; dropped. First few: {sorted(transcriptomic_only)[:10]}"
        )

    fused = clinical_encoded.loc[common_ids].join(hybrid_scores.loc[common_ids], how="inner")
    fused.index.name = "sample_id"

    print(
        f"[FUSION][{cohort_label}] Fused {fused.shape[0]} samples x {fused.shape[1]} features "
        f"({clinical_encoded.shape[1]} clinical + {hybrid_scores.shape[1]} transcriptomic)."
    )

    return fused

CLINICAL_PATH      = "data/raw/combined_metadata.csv"
TRAIN_SCORES_PATH  = "outputs/scores/gse_train_scores.csv"
TEST_SCORES_PATH   = "outputs/scores/gse_test_scores.csv"
TRAIN_LE_PATH       = "outputs/scores/gse_train_leading_edge.csv"
TEST_LE_PATH        = "outputs/scores/gse_test_leading_edge.csv"
TRAIN_TF_PATH       = "outputs/scores/gse_train_tf_scores.csv"
TEST_TF_PATH        = "outputs/scores/gse_test_tf_scores.csv"
FUSED_OUTPUT_DIR    = "outputs/python"


def fuse_all_cohorts() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fuse train and test cohorts only. No ICGA. Verifies the two fused
    DataFrames have identical column sets -- if they don't (e.g. a
    categorical level that only appears in one partition despite the fixed
    schema, or a dropped-sample edge case), aligns both to the column
    intersection and warns."""

    train_fused = fuse_features(CLINICAL_PATH, TRAIN_SCORES_PATH, TRAIN_LE_PATH, TRAIN_TF_PATH,
                                 cohort_label="train")
    test_fused  = fuse_features(CLINICAL_PATH, TEST_SCORES_PATH, TEST_LE_PATH, TEST_TF_PATH,
                                 cohort_label="test")

    train_cols = set(train_fused.columns)
    test_cols  = set(test_fused.columns)

    if train_cols != test_cols:
        only_train = train_cols - test_cols
        only_test  = test_cols - train_cols
        warnings.warn(
            f"[FUSION] Train/test column sets differ -- aligning to the intersection. "
            f"{len(only_train)} column(s) only in train ({sorted(only_train)[:10]}), "
            f"{len(only_test)} column(s) only in test ({sorted(only_test)[:10]})."
        )

    common_cols = sorted(train_cols & test_cols)
    train_fused = train_fused[common_cols]
    test_fused  = test_fused[common_cols]

    out_dir = Path(FUSED_OUTPUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    train_out = out_dir / "fused_gse_train.csv"
    test_out  = out_dir / "fused_gse_test.csv"

    train_fused.to_csv(train_out)
    test_fused.to_csv(test_out)

    print(
        f"[FUSION] Saved {train_out} ({train_fused.shape[0]}x{train_fused.shape[1]}) and "
        f"{test_out} ({test_fused.shape[0]}x{test_fused.shape[1]})."
    )

    return train_fused, test_fused