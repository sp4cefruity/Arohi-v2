"""
AROHI v2 -- run_pipeline.py

Orchestrates Phase 4 (feature fusion + selection) and Phase 5 (model
training) end to end. The CONFIG block below is the ONLY part you should
need to edit between runs -- everything else in this file is fixed
plumbing.

Ablation for cohort-agnostic modeling:
    python run_pipeline.py --ablate
runs three variants sequentially (all features / gse_id dropped /
gse_id dropped + top-150 MI features) and prints a comparison of the
honest CV AUCs. Single variants can be run with:
    python run_pipeline.py --drop-gse-id [--k 150] [--n-trials 60]
"""

from __future__ import annotations

import argparse
import shutil
import warnings
from pathlib import Path

import pandas as pd

from evaluation import compute_metrics_dual, load_model_bundle, predict_cohort_dual
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

# AUC-optimization knobs (see --ablate):
#   drop_gse_id      -- remove the 31 gse_id_* study one-hots so the models
#                       can't memorize study-level pCR rates.
#   feature_reduction_k -- keep only the top-k mutual-information features
#                       for BOTH the ElasticNet selection and XGBoost.
#   tune_cv_folds / n_trials -- Optuna search cost.
#   report_cv_folds   -- folds for the honest reported metrics (seed 123,
#                       distinct from the tuning folds).
DROP_GSE_ID        = False
FEATURE_REDUCTION_K = None
N_TRIALS            = 60
TUNE_CV_FOLDS       = 5
REPORT_CV_FOLDS     = 10
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


def run_once(drop_gse_id: bool, feature_reduction_k: int | None,
             n_trials: int, tune_cv_folds: int, report_cv_folds: int) -> dict:
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
        cv_folds=SELECTION_CV_FOLDS, random_state=RANDOM_STATE, output_dir=OUTPUT_DIR,
        drop_gse_id=drop_gse_id, feature_reduction_k=feature_reduction_k,
    )

    print("[PIPELINE] Model training (GLM + XGBoost + stacking)...")
    training_result = run_model_training_pipeline(
        selection_result,
        cv_folds=TRAINING_CV_FOLDS, output_dir=OUTPUT_DIR,
        tune_cv_folds=tune_cv_folds, report_cv_folds=report_cv_folds,
        n_trials=n_trials, drop_gse_id=drop_gse_id, feature_reduction_k=feature_reduction_k,
    )

    print("[PIPELINE] Done.")

    return {
        "fused_train":      fused_train,
        "fused_test":       fused_test,
        "selection_result": selection_result,
        "training_result":  training_result,
        "config":           {
            "drop_gse_id": drop_gse_id,
            "feature_reduction_k": feature_reduction_k,
            "n_trials": n_trials,
            "tune_cv_folds": tune_cv_folds,
            "report_cv_folds": report_cv_folds,
        },
    }


def evaluate_on_test(bundle_path: str, output_dir: str) -> dict:
    """Lightweight held-out test evaluation (no SHAP / plots): load the
    bundle, score the fused test set, return per-model test AUC."""
    bundle = load_model_bundle(bundle_path)
    fused_test = pd.read_csv(f"{output_dir}/fused_gse_test.csv", index_col=0)
    fused_test.index = fused_test.index.astype(str)
    fused_test.index.name = "sample_id"
    labels = load_labels(LABELS_TRAIN_PATH)

    preds = predict_cohort_dual(bundle, fused_test, "test")
    common = preds.index.intersection(labels.index)
    metrics = compute_metrics_dual(labels.loc[common], preds.loc[common], "test")
    return {
        "xgb_test_auc":  metrics["xgb"]["auc"],
        "ens_test_auc":  metrics["ensemble"]["auc"],
        "ens_test_auprc": metrics["ensemble"]["auprc"],
    }


def print_ablation_summary(results: list[dict]):
    print("\n==================== ABLATION SUMMARY ====================")
    print(f"{'variant':<26}{'glm AUC':>8}{'xgb AUC':>8}{'stack AUC':>9}"
          f"{'xgb test AUC':>13}{'ens test AUC':>13}")
    for r in results:
        meta = r["training_result"]["training_metadata"]
        variant = "gse_id kept"
        if meta["drop_gse_id"]:
            variant = "gse_id dropped"
        if meta["feature_reduction_k"]:
            variant = f"gse_id dropped, top-{meta['feature_reduction_k']} MI"
        test = r["test_metrics"]
        print(
            f"{variant:<26}"
            f"{meta['glm_cv_auc_mean']:>8.3f}"
            f"{meta['xgb_cv_auc_mean']:>8.3f}"
            f"{meta['stack_cv_auc_mean']:>9.3f}"
            f"{test['xgb_test_auc']:>13.3f}"
            f"{test['ens_test_auc']:>13.3f}"
        )
    print("=" * 74)


def main() -> dict:
    parser = argparse.ArgumentParser(description="AROHI v2 Phase 4/5 pipeline")
    parser.add_argument("--drop-gse-id", action="store_true",
                        help="Drop the 31 gse_id study one-hot columns before modeling.")
    parser.add_argument("--k", type=int, default=None,
                        help="Keep only the top-k mutual-information features.")
    parser.add_argument("--n-trials", type=int, default=N_TRIALS,
                        help="Optuna trials for XGBoost tuning.")
    parser.add_argument("--tune-folds", type=int, default=TUNE_CV_FOLDS)
    parser.add_argument("--report-folds", type=int, default=REPORT_CV_FOLDS)
    parser.add_argument("--ablate", action="store_true",
                        help="Run all three variants and print a comparison.")
    args = parser.parse_args()

    if args.ablate:
        variants = [
            {"drop_gse_id": False, "feature_reduction_k": None},
            {"drop_gse_id": True,  "feature_reduction_k": None},
            {"drop_gse_id": True,  "feature_reduction_k": 150},
        ]
        results = []
        for i, v in enumerate(variants, start=1):
            print(f"\n########## VARIANT: drop_gse_id={v['drop_gse_id']}, "
                  f"feature_reduction_k={v['feature_reduction_k']} ##########")
            res = run_once(
                drop_gse_id=v["drop_gse_id"], feature_reduction_k=v["feature_reduction_k"],
                n_trials=args.n_trials, tune_cv_folds=args.tune_folds,
                report_cv_folds=args.report_folds,
            )
            bundle_src = Path(OUTPUT_DIR) / "model_bundle.joblib"
            bundle_dst = Path(OUTPUT_DIR) / f"model_bundle_v{i}.joblib"
            shutil.copy(bundle_src, bundle_dst)
            res["test_metrics"] = evaluate_on_test(str(bundle_dst), OUTPUT_DIR)
            results.append(res)
        print_ablation_summary(results)
        return {"results": results}

    return run_once(
        drop_gse_id=args.drop_gse_id,
        feature_reduction_k=args.k,
        n_trials=args.n_trials,
        tune_cv_folds=args.tune_folds,
        report_cv_folds=args.report_folds,
    )


if __name__ == "__main__":
    main()
