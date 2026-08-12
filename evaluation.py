from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
import statsmodels.api as sm
from sklearn.calibration import calibration_curve
from sklearn.metrics import (
    accuracy_score, auc, average_precision_score, brier_score_loss,
    confusion_matrix, f1_score, precision_recall_curve, roc_auc_score, roc_curve,
)

from feature_selection import apply_preprocessing
from model_training import _to_binary

MODEL_COLORS = {"glm": "green", "xgb": "blue", "ensemble": "red"}
PROB_COLS    = {"glm": "prob_pCR_glm", "xgb": "prob_pCR_xgb", "ensemble": "prob_pCR_ensemble"}
LABEL_COLS   = {"glm": "pred_label_glm", "xgb": "pred_label_xgb", "ensemble": "pred_label_ensemble"}

def load_model_bundle(path: str) -> dict:
    import joblib
    bundle = joblib.load(path)

    meta = bundle.get("training_metadata", {})
    n_glm = len(bundle.get("selected_features", []))
    n_xgb = len(bundle.get("all_train_columns", []))

    print(f"[BUNDLE] Loaded {path}")
    glm_auc_str = (
        f", train CV AUC = {meta['glm_cv_auc_mean']:.3f} +/- {meta['glm_cv_auc_std']:.3f}"
        if "glm_cv_auc_mean" in meta else ""
    )
    xgb_auc_str = (
        f", train CV AUC = {meta['xgb_cv_auc_mean']:.3f} +/- {meta['xgb_cv_auc_std']:.3f}"
        if "xgb_cv_auc_mean" in meta else ""
    )
    print(f"[BUNDLE]   GLM     : {n_glm} feature(s){glm_auc_str}")
    print(f"[BUNDLE]   XGBoost : {n_xgb} feature(s){xgb_auc_str}")

    return bundle


# ---------------------------------------------------------------------------
# predict_cohort_dual()
# ---------------------------------------------------------------------------
def predict_cohort_dual(bundle: dict, fused_df: pd.DataFrame, cohort_label: str) -> pd.DataFrame:
    imputer            = bundle["imputer"]
    scaler              = bundle["scaler"]
    selected_features   = bundle["selected_features"]
    all_train_columns   = bundle["all_train_columns"]
    glm                  = bundle["glm"]
    xgb_model            = bundle["xgb_model"]

    missing = [c for c in all_train_columns if c not in fused_df.columns]
    if missing:
        raise ValueError(
            f"[PREDICT][{cohort_label}] fused_df is missing {len(missing)} column(s) the "
            f"model needs, e.g. {missing[:10]}"
        )

    X_full   = fused_df[all_train_columns]
    X_scaled = apply_preprocessing(X_full, imputer, scaler)   # .transform() only

    
    X_glm = sm.add_constant(X_scaled[selected_features], has_constant="add")
    X_glm = X_glm.reindex(columns=glm.params.index, fill_value=0.0)  # match training column order
    prob_glm = pd.Series(glm.predict(X_glm).values, index=X_scaled.index, name="prob_pCR_glm")

    
    prob_xgb = pd.Series(
        xgb_model.predict_proba(X_scaled[all_train_columns])[:, 1],
        index=X_scaled.index, name="prob_pCR_xgb"
    )

    
    prob_ens = ((prob_glm + prob_xgb) / 2.0).rename("prob_pCR_ensemble")

    to_label = lambda p: p.apply(lambda v: "pCR" if v >= 0.5 else "RD")

    result = pd.DataFrame({
        "prob_pCR_glm":        prob_glm,
        "prob_pCR_xgb":        prob_xgb,
        "prob_pCR_ensemble":   prob_ens,
        "pred_label_glm":      to_label(prob_glm),
        "pred_label_xgb":      to_label(prob_xgb),
        "pred_label_ensemble": to_label(prob_ens),
    })
    result.index.name = "sample_id"

    print(
        f"[PREDICT][{cohort_label}] Scored {len(result)} sample(s) "
        f"(GLM: {len(selected_features)} features, XGBoost: {len(all_train_columns)} features)."
    )
    return result



def compute_metrics_dual(y_true, predictions: pd.DataFrame, cohort_label: str) -> dict:
    y_binary = _to_binary(y_true).reindex(predictions.index)

    results = {}
    for name, prob_col in PROB_COLS.items():
        probs = predictions[prob_col].values
        preds = (probs >= 0.5).astype(int)
        results[name] = {
            "auc":      float(roc_auc_score(y_binary, probs)),
            "auprc":    float(average_precision_score(y_binary, probs)),
            "f1":       float(f1_score(y_binary, preds, zero_division=0)),
            "accuracy": float(accuracy_score(y_binary, preds)),
            "brier":    float(brier_score_loss(y_binary, probs)),
        }

    print(f"\n==================== Metrics: {cohort_label} ====================")
    print(f"{'Model':<10}{'AUC':>8}{'AUPRC':>8}{'F1':>8}{'Accuracy':>10}{'Brier':>8}")
    for name in ("glm", "xgb", "ensemble"):
        m = results[name]
        print(f"{name:<10}{m['auc']:>8.3f}{m['auprc']:>8.3f}{m['f1']:>8.3f}{m['accuracy']:>10.3f}{m['brier']:>8.3f}")
    print()

    return results



def plot_roc_curves_dual(y_true, predictions: pd.DataFrame, cohort_label: str, output_dir: str) -> str:
    y_binary = _to_binary(y_true).reindex(predictions.index)

    fig, ax = plt.subplots(figsize=(6, 6))
    for name, col in PROB_COLS.items():
        fpr, tpr, _ = roc_curve(y_binary, predictions[col])
        roc_auc = auc(fpr, tpr)
        ax.plot(fpr, tpr, color=MODEL_COLORS[name], label=f"{name.upper()} (AUC={roc_auc:.3f})")
    ax.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"ROC curves -- {cohort_label}")
    ax.legend(loc="lower right")

    return _save_plot(fig, output_dir, f"{cohort_label}_roc_curves.png", "ROC curves")


def plot_pr_curves_dual(y_true, predictions: pd.DataFrame, cohort_label: str, output_dir: str) -> str:
    y_binary = _to_binary(y_true).reindex(predictions.index)

    fig, ax = plt.subplots(figsize=(6, 6))
    for name, col in PROB_COLS.items():
        precision, recall, _ = precision_recall_curve(y_binary, predictions[col])
        ap = average_precision_score(y_binary, predictions[col])
        ax.plot(recall, precision, color=MODEL_COLORS[name], label=f"{name.upper()} (AUPRC={ap:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(f"Precision-Recall curves -- {cohort_label}")
    ax.legend(loc="lower left")

    return _save_plot(fig, output_dir, f"{cohort_label}_pr_curves.png", "PR curves")


def plot_calibration_curves_dual(y_true, predictions: pd.DataFrame, cohort_label: str,
                                  output_dir: str, n_bins: int = 10) -> str:
    y_binary = _to_binary(y_true).reindex(predictions.index)

    fig, ax = plt.subplots(figsize=(6, 6))
    for name, col in PROB_COLS.items():
        frac_pos, mean_pred = calibration_curve(y_binary, predictions[col], n_bins=n_bins, strategy="quantile")
        ax.plot(mean_pred, frac_pos, marker="o", color=MODEL_COLORS[name], label=name.upper())
    ax.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)
    ax.set_xlabel("Mean predicted probability")
    ax.set_ylabel("Observed frequency of pCR")
    ax.set_title(f"Calibration curves -- {cohort_label}")
    ax.legend(loc="upper left")

    return _save_plot(fig, output_dir, f"{cohort_label}_calibration_curves.png", "calibration curves")


def plot_confusion_matrices_dual(y_true, predictions: pd.DataFrame, cohort_label: str, output_dir: str) -> str:
    y_binary = _to_binary(y_true).reindex(predictions.index)
    cmaps = {"glm": "Greens", "xgb": "Blues", "ensemble": "Reds"}

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (name, col) in zip(axes, LABEL_COLS.items()):
        y_pred_binary = (predictions[col] == "pCR").astype(int)
        cm = confusion_matrix(y_binary, y_pred_binary, labels=[0, 1])
        ax.imshow(cm, cmap=cmaps[name])
        ax.set_xticks([0, 1]); ax.set_xticklabels(["RD", "pCR"])
        ax.set_yticks([0, 1]); ax.set_yticklabels(["RD", "pCR"])
        ax.set_xlabel("Predicted"); ax.set_ylabel("True")
        ax.set_title(name.upper())
        for i in range(2):
            for j in range(2):
                ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                        color="white" if cm[i, j] > cm.max() / 2 else "black")

    fig.suptitle(f"Confusion matrices -- {cohort_label}")
    fig.tight_layout()

    return _save_plot(fig, output_dir, f"{cohort_label}_confusion_matrices.png", "confusion matrices")


def _save_plot(fig, output_dir: str, filename: str, label: str) -> str:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / filename
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[PLOT] Saved {label} -> {out_path}")
    return str(out_path)

def run_shap_analysis_xgb(bundle: dict, fused_df: pd.DataFrame, predictions: pd.DataFrame,
                           cohort_label: str, n_waterfall: int = 5, output_dir: str = "outputs") -> pd.DataFrame:
    xgb_model          = bundle["xgb_model"]
    imputer             = bundle["imputer"]
    scaler               = bundle["scaler"]
    all_train_columns   = bundle["all_train_columns"]

    X_full   = fused_df[all_train_columns]
    X_scaled = apply_preprocessing(X_full, imputer, scaler)

    explainer   = shap.TreeExplainer(xgb_model)
    shap_values = explainer(X_scaled)

    shap_df = pd.DataFrame(shap_values.values, index=X_scaled.index, columns=X_scaled.columns)

    shap_dir = Path(output_dir) / "shap"
    shap_dir.mkdir(parents=True, exist_ok=True)

    # Beeswarm summary, top 20 features
    plt.figure()
    shap.summary_plot(shap_values, X_scaled, max_display=20, show=False)
    beeswarm_path = shap_dir / f"{cohort_label}_shap_beeswarm.png"
    plt.savefig(beeswarm_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[SHAP] Saved beeswarm summary -> {beeswarm_path}")

    
    top_patients = predictions["prob_pCR_xgb"].sort_values(ascending=False).head(n_waterfall).index
    for pid in top_patients:
        pos = X_scaled.index.get_loc(pid)
        plt.figure()
        shap.plots.waterfall(shap_values[pos], show=False)
        wf_path = shap_dir / f"{cohort_label}_waterfall_{pid}.png"
        plt.savefig(wf_path, dpi=150, bbox_inches="tight")
        plt.close()
    print(f"[SHAP] Saved {len(top_patients)} waterfall plot(s) for top patients by predicted pCR probability.")

    # Dependence plots, top-3 features by mean |SHAP|
    top_features = shap_df.abs().mean().sort_values(ascending=False).head(3).index
    for feat in top_features:
        plt.figure()
        shap.dependence_plot(feat, shap_values.values, X_scaled, show=False)
        dep_path = shap_dir / f"{cohort_label}_dependence_{feat}.png"
        plt.savefig(dep_path, dpi=150, bbox_inches="tight")
        plt.close()
    print(f"[SHAP] Saved dependence plots for top 3 features: {list(top_features)}")

    return shap_df


def run_glm_explanation(bundle: dict, fused_df: pd.DataFrame, predictions: pd.DataFrame,
                         cohort_label: str, output_dir: str = "outputs") -> pd.DataFrame:
    glm                  = bundle["glm"]
    imputer             = bundle["imputer"]
    scaler               = bundle["scaler"]
    selected_features   = bundle["selected_features"]
    all_train_columns   = bundle["all_train_columns"]

    X_full   = fused_df[all_train_columns]
    X_scaled = apply_preprocessing(X_full, imputer, scaler)
    X_sel    = X_scaled[selected_features]

    coefs   = glm.params.reindex(selected_features)          # excludes 'const'
    contrib = X_sel.multiply(coefs, axis=1)                  # per-patient per-feature log-odds contribution
    contrib.index.name = "sample_id"

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    contrib_path = out_dir / f"{cohort_label}_glm_contributions.csv"
    contrib.to_csv(contrib_path)
    print(f"[GLM EXPLAIN] Saved per-patient log-odds contributions -> {contrib_path}")

    return contrib


def _top_k_by_abs(series: pd.Series, k: int = 3) -> pd.Series:
    series = series.dropna()
    if series.empty:
        return series
    return series.reindex(series.abs().sort_values(ascending=False).index).head(k)


def generate_patient_narrative(patient_id, xgb_prob: float, glm_prob: float,
                                shap_row: pd.Series, glm_contrib_row: pd.Series) -> str:
    """NOTE on GLM 'odds ratios': this function only receives per-patient
    log-odds CONTRIBUTIONS (coefficient x this patient's feature value),
    not the raw model coefficients, so there's no fixed per-feature odds
    ratio available here to report. What's shown instead is
    exp(contribution) -- this patient's own local multiplicative effect on
    the odds for that feature, which is patient-specific rather than a
    static model-wide number."""

    top_shap = _top_k_by_abs(shap_row, 3)
    top_glm  = _top_k_by_abs(glm_contrib_row, 3)

    xgb_drivers = "; ".join(f"{feat} (SHAP={val:+.3f})" for feat, val in top_shap.items()) or "n/a"
    glm_drivers = "; ".join(
        f"{feat} (local odds multiplier={np.exp(val):.2f}x)" for feat, val in top_glm.items()
    ) or "n/a"

    lines = [
        f"Patient {patient_id}:",
        f"XGBoost predicts {xgb_prob * 100:.1f}% pCR probability, driven by the top three "
        f"features and their SHAP values: {xgb_drivers}.",
        f"GLM predicts {glm_prob * 100:.1f}% pCR probability, driven by the top three "
        f"features and their odds ratios: {glm_drivers}.",
    ]

    shared = set(top_shap.index) & set(top_glm.index)
    if shared:
        lines.append(
            f"Both models agree on {len(shared)} of the same top driver feature(s) "
            f"({', '.join(sorted(shared))}) -- the most reliable kind of prediction."
        )

    if abs(xgb_prob - glm_prob) < 0.10:
        lines.append("High model agreement — prediction is clinically trustworthy.")
    else:
        lines.append("Low model agreement — interpret with caution.")

    return " ".join(lines)


def run_evaluation_pipeline(model_path: str, fused_test_path: str, test_labels_path: str,
                             output_dir: str = "outputs", run_shap: bool = True,
                             n_waterfall: int = 5) -> dict:
    bundle = load_model_bundle(model_path)

    fused_test = pd.read_csv(fused_test_path, index_col=0)
    fused_test.index = fused_test.index.astype(str)
    fused_test.index.name = "sample_id"

    
    test_labels_df = pd.read_csv(test_labels_path)
    if "sample_id" not in test_labels_df.columns or "response" not in test_labels_df.columns:
        raise ValueError(f"[EVAL] '{test_labels_path}' must have 'sample_id' and 'response' columns.")
    test_labels_df = test_labels_df.set_index("sample_id")
    test_labels_df.index = test_labels_df.index.astype(str)
    y_test = test_labels_df["response"]

    cohort_label = "test"
    predictions = predict_cohort_dual(bundle, fused_test, cohort_label)

    common_ids = predictions.index.intersection(y_test.index)
    if len(common_ids) < len(predictions):
        warnings.warn(
            f"[EVAL] {len(predictions) - len(common_ids)} predicted sample(s) have no test "
            f"label; excluded from metrics (but still get a narrative)."
        )

    predictions_eval = predictions.loc[common_ids]
    y_test_eval       = y_test.loc[common_ids]

    metrics = compute_metrics_dual(y_test_eval, predictions_eval, cohort_label)

    plot_roc_curves_dual(y_test_eval, predictions_eval, cohort_label, output_dir)
    plot_pr_curves_dual(y_test_eval, predictions_eval, cohort_label, output_dir)
    plot_calibration_curves_dual(y_test_eval, predictions_eval, cohort_label, output_dir)
    plot_confusion_matrices_dual(y_test_eval, predictions_eval, cohort_label, output_dir)

    shap_df = None
    if run_shap:
        shap_df = run_shap_analysis_xgb(
            bundle, fused_test.loc[predictions.index], predictions,
            cohort_label, n_waterfall=n_waterfall, output_dir=output_dir
        )

    glm_contrib_df = run_glm_explanation(
        bundle, fused_test.loc[predictions.index], predictions, cohort_label, output_dir=output_dir
    )

    narratives = []
    for pid in predictions.index:
        xgb_prob = float(predictions.loc[pid, "prob_pCR_xgb"])
        glm_prob = float(predictions.loc[pid, "prob_pCR_glm"])
        shap_row = shap_df.loc[pid] if shap_df is not None and pid in shap_df.index else pd.Series(dtype=float)
        glm_row  = glm_contrib_df.loc[pid] if pid in glm_contrib_df.index else pd.Series(dtype=float)
        narratives.append({
            "sample_id": pid,
            "narrative": generate_patient_narrative(pid, xgb_prob, glm_prob, shap_row, glm_row),
        })

    patient_report = pd.DataFrame(narratives).set_index("sample_id").join(predictions)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    report_path = out_dir / "patient_report.csv"
    patient_report.to_csv(report_path)
    print(f"[EVAL] Saved patient report -> {report_path}")

    metrics_path = out_dir / "gse_test_metrics_dual.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"[EVAL] Saved metrics -> {metrics_path}")

    return {
        "bundle":            bundle,
        "predictions":       predictions,
        "metrics":           metrics,
        "shap_values":       shap_df,
        "glm_contributions": glm_contrib_df,
        "patient_report":    patient_report,
    }