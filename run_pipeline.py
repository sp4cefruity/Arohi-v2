"""
AROHI v2 -- run_pipeline.py

Orchestrates Phase 4 (feature fusion + selection) and Phase 5 (model
training) end to end. The CONFIG block below is the ONLY part you should
need to edit between runs -- everything else in this file is fixed
plumbing.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd

from feature_fusion import fuse_features
from feature_selection import run_feature_selection_pipeline
from model_training import run_model_training_pipeline

# =============================================================================
# CONFIG -- edit these between runs. Nothing below this block should need
# to change.
# =============================================================================
GSE_TRAIN_CLINICAL = "data/raw/combined_metadata.csv"
GSE_TRAIN_SCORES   = "outputs/scores/gse_train_ssgsea.csv"
GSE_TRAIN_LE       = "outputs/scores/gse_train_leading_edge.csv"
GSE_TRAIN_TF       = "outputs/scores/gse_train_tf_scores.csv"

GSE_TEST_CLINICAL  = "data/raw/combined_metadata.csv"
GSE_TEST_SCORES    = "outputs/scores/gse_test_ssgsea.csv"
GSE_TEST_LE        = "outputs/scores/gse_test_leading_edge.csv"
GSE_TEST_TF        = "outputs/scores/gse_test_tf_scores.csv"

LABELS_TRAIN_PATH  = "data/raw/combined_labels.csv"

OUTPUT_DIR          = "outputs/python"

SELECTION_CV_FOLDS = 10
TRAINING_CV_FOLDS   = 10
RANDOM_STATE         = 42
# =============================================================================


def load_labels(labels_path: str) -> pd.Series:
    """data/raw/combined_labels.csv -> a sample_id-indexed Series of
    'pCR'/'RD' response labels."""
    df = pd.read_csv(labels_path)
    if "sample_id" not in df.columns or "response" not in df.columns:
        raise ValueError(f"[LABELS] '{labels_path}' must have 'sample_id' and 'response' columns.")
    df = df.set_index("sample_id")
    df.index = df.index.astype(str)
    return df["response"]


def main() -> dict:
    print("[PIPELINE] === AROHI v2 -- Phase 4/5 (Python) ===")

    print("[PIPELINE] Fusing train cohort...")
    fused_train = fuse_features(
        GSE_TRAIN_CLINICAL, GSE_TRAIN_SCORES, GSE_TRAIN_LE, GSE_TRAIN_TF, cohort_label="train"
    )

    print("[PIPELINE] Fusing test cohort...")
    fused_test = fuse_features(
        GSE_TEST_CLINICAL, GSE_TEST_SCORES, GSE_TEST_LE, GSE_TEST_TF, cohort_label="test"
    )

    # Align train/test to a common column set (same guard as
    # feature_fusion.fuse_all_cohorts(), inlined here since these paths are
    # config-driven rather than that function's hardcoded defaults).
    train_cols = set(fused_train.columns)
    test_cols  = set(fused_test.columns)
    if train_cols != test_cols:
        only_train = train_cols - test_cols
        only_test  = test_cols - train_cols
        warnings.warn(
            f"[PIPELINE] Train/test column sets differ -- aligning to the intersection. "
            f"{len(only_train)} column(s) only in train, {len(only_test)} column(s) only in test."
        )
    common_cols = sorted(train_cols & test_cols)
    fused_train = fused_train[common_cols]
    fused_test  = fused_test[common_cols]

    out_dir = Path(OUTPUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    fused_train.to_csv(out_dir / "fused_gse_train.csv")
    fused_test.to_csv(out_dir / "fused_gse_test.csv")
    print(f"[PIPELINE] Saved fused_gse_train.csv ({fused_train.shape}) and fused_gse_test.csv ({fused_test.shape}).")

    labels_train = load_labels(LABELS_TRAIN_PATH)

    print("[PIPELINE] Feature selection (ElasticNetCV)...")
    selection_result = run_feature_selection_pipeline(
        fused_train, labels_train,
        cv_folds=SELECTION_CV_FOLDS, random_state=RANDOM_STATE, output_dir=OUTPUT_DIR
    )

    print("[PIPELINE] Model training (GLM + XGBoost)...")
    training_result = run_model_training_pipeline(
        selection_result, cv_folds=TRAINING_CV_FOLDS, output_dir=OUTPUT_DIR
    )

    print("[PIPELINE] Done.")

    return {
        "fused_train":      fused_train,
        "fused_test":       fused_test,
        "selection_result": selection_result,
        "training_result":  training_result,
    }


if __name__ == "__main__":
    main()