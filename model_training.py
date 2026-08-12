from __future__ import annotations

import warnings
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.base import clone as sk_clone
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from xgboost import XGBClassifier

FEATURE_SELECTION_SEED = 42   # must stay distinct from XGBOOST_TUNING_SEED
XGBOOST_TUNING_SEED     = 7


def compute_scale_pos_weight(y_train) -> float:
    """n_RD / n_pCR, for XGBoost's scale_pos_weight. Accepts either string
    labels ('pCR'/'RD') or a binary encoding (1=pCR, 0=RD)."""
    s = pd.Series(y_train)

    if s.dtype == object or str(s.dtype).startswith("string") or str(s.dtype) == "category":
        s_str = s.astype(str).str.strip()
        n_pcr = int((s_str == "pCR").sum())
        n_rd  = int((s_str == "RD").sum())
    else:
        n_pcr = int((s == 1).sum())
        n_rd  = int((s == 0).sum())

    if n_pcr == 0:
        raise ValueError("[SCALE_POS_WEIGHT] No pCR samples found; cannot compute scale_pos_weight.")

    weight  = n_rd / n_pcr
    pct_pcr = 100 * n_pcr / (n_pcr + n_rd)
    print(
        f"[CLASS BALANCE] n_pCR={n_pcr}, n_RD={n_rd} ({pct_pcr:.1f}% pCR) "
        f"-> scale_pos_weight={weight:.4f}"
    )
    return weight


def train_glm(X_train: pd.DataFrame, y_binary, selected_features: list[str],
              output_dir: str = "outputs"):
    """sm.Logit on the ElasticNet-selected features + intercept. Convergence
    warnings are caught and logged, not allowed to crash the run."""

    X_design = sm.add_constant(X_train[selected_features], has_constant="add")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        try:
            model = sm.Logit(y_binary, X_design).fit(disp=0, maxiter=500)
        except Exception as e:
            raise RuntimeError(f"[GLM] sm.Logit failed to fit: {e}") from e

    for w in caught:
        print(f"[GLM][WARN] {w.message}")

    converged = getattr(model, "mle_retvals", {}).get("converged", None)
    if converged is False:
        print("[GLM][WARN] Did not converge within maxiter=500; coefficients may be unreliable.")

    coef_table = pd.DataFrame({
        "feature":     X_design.columns,
        "coefficient": model.params.values,
        "p_value":     model.pvalues.values,
    })
    coef_table["odds_ratio"] = np.exp(coef_table["coefficient"])

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "glm_coefficients.csv"
    coef_table.to_csv(out_path, index=False)

    print(
        f"[GLM] Fit logistic regression on {len(selected_features)} selected feature(s) "
        f"+ intercept -> {out_path}"
    )

    return model


def _to_binary(y) -> pd.Series:
    """pCR -> 1, RD -> 0. Passes through an already-binary series unchanged."""
    s = pd.Series(y)
    if s.dtype == object or str(s.dtype).startswith("string") or str(s.dtype) == "category":
        return (s.astype(str).str.strip() == "pCR").astype(int)
    return s.astype(int)


def tune_xgboost(X_train: pd.DataFrame, y_binary, scale_pos_weight: float, cv_folds: int = 10):
    """GridSearchCV over the specified grid, refit on all training data with
    the best params. Uses XGBOOST_TUNING_SEED (7) -- deliberately different
    from feature selection's seed (42)."""

    y = np.asarray(pd.Series(y_binary))
    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=XGBOOST_TUNING_SEED)

    param_grid = {
        "max_depth":        [2, 3, 4],
        "learning_rate":    [0.01, 0.05, 0.1],
        "n_estimators":     [100, 200, 300],
        "subsample":        [0.7, 0.9],
        "colsample_bytree": [0.7, 0.9],
        "min_child_weight": [1, 3, 5],
    }

    base_model = XGBClassifier(
        objective="binary:logistic",
        eval_metric="auc",
        scale_pos_weight=scale_pos_weight,
        n_jobs=1,
        random_state=XGBOOST_TUNING_SEED,
    )

    grid_search = GridSearchCV(
        estimator=base_model,
        param_grid=param_grid,
        scoring="roc_auc",
        cv=cv,
        refit=True,
        n_jobs=-1,
    )
    grid_search.fit(X_train, y)

    print(f"[XGB TUNE] Best CV AUC: {grid_search.best_score_:.4f}")
    print(f"[XGB TUNE] Best params: {grid_search.best_params_}")

    return grid_search.best_estimator_


def cross_validated_metrics(model, X_train: pd.DataFrame, y_binary, cv_folds: int = 10,
                             model_name: str = "model", output_dir: str = "outputs") -> pd.DataFrame:
    """Out-of-fold AUC / AUPRC / F1 / accuracy. Refits a fresh clone of
    `model` inside each fold (sklearn estimators via sklearn.base.clone;
    statsmodels Logit results via a fresh sm.Logit call using the same
    exogenous columns) so metrics aren't leaked from the already-fit model
    passed in. Uses the same random_state=7 folds as XGBoost tuning."""

    y = np.asarray(pd.Series(y_binary))
    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=XGBOOST_TUNING_SEED)

    is_sklearn = hasattr(model, "get_params") and hasattr(model, "fit") and hasattr(model, "predict_proba")
    if not is_sklearn:
        exog_names = [c for c in model.model.exog_names if c != "const"]

    rows = []
    for fold_idx, (train_idx, val_idx) in enumerate(cv.split(X_train, y), start=1):
        X_fold_train, X_fold_val = X_train.iloc[train_idx], X_train.iloc[val_idx]
        y_fold_train, y_fold_val = y[train_idx], y[val_idx]

        if is_sklearn:
            fold_model = sk_clone(model)
            fold_model.fit(X_fold_train, y_fold_train)
            y_pred_proba = fold_model.predict_proba(X_fold_val)[:, 1]
        else:
            design_train = sm.add_constant(X_fold_train[exog_names], has_constant="add")
            design_val   = sm.add_constant(X_fold_val[exog_names], has_constant="add")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", ConvergenceWarning)
                fold_model = sm.Logit(y_fold_train, design_train).fit(disp=0, maxiter=500)
            y_pred_proba = fold_model.predict(design_val)

        y_pred_label = (y_pred_proba >= 0.5).astype(int)

        rows.append({
            "fold":     fold_idx,
            "auc":      roc_auc_score(y_fold_val, y_pred_proba),
            "auprc":    average_precision_score(y_fold_val, y_pred_proba),
            "f1":       f1_score(y_fold_val, y_pred_label, zero_division=0),
            "accuracy": accuracy_score(y_fold_val, y_pred_label),
        })

    results = pd.DataFrame(rows)
    summary = results[["auc", "auprc", "f1", "accuracy"]].agg(["mean", "std"])
    print(
        f"[CV METRICS][{model_name}] "
        f"AUC={summary.loc['mean','auc']:.3f}±{summary.loc['std','auc']:.3f}  "
        f"AUPRC={summary.loc['mean','auprc']:.3f}±{summary.loc['std','auprc']:.3f}  "
        f"F1={summary.loc['mean','f1']:.3f}±{summary.loc['std','f1']:.3f}  "
        f"Accuracy={summary.loc['mean','accuracy']:.3f}±{summary.loc['std','accuracy']:.3f}"
    )

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"cv_results_{model_name}.csv"
    results.to_csv(out_path, index=False)
    print(f"[CV METRICS][{model_name}] Saved -> {out_path}")

    return results


def save_model_bundle(glm, xgb_model, imputer, scaler, selected_features, all_train_columns,
                       training_metadata: dict, output_dir: str = "outputs/python") -> str:
    """Single joblib dump containing everything needed to reproduce
    predictions with no refitting."""
    bundle = {
        "glm":                glm,
        "xgb_model":          xgb_model,
        "imputer":            imputer,
        "scaler":             scaler,
        "selected_features":  selected_features,
        "all_train_columns":  all_train_columns,
        "training_metadata":  training_metadata,
    }

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "model_bundle.joblib"
    joblib.dump(bundle, out_path)

    print(f"[BUNDLE] Saved model bundle -> {out_path}")
    return str(out_path)


def plot_feature_importance_dual(glm, xgb_model, selected_features, all_columns,
                                  output_dir: str = "outputs/python") -> str:
    """Side-by-side: GLM odds ratios (left) vs XGBoost gain importance
    (right), top 20 features each, ranked by |log-odds| and gain
    respectively."""

    glm_params = glm.params.drop(labels=["const"], errors="ignore")
    order      = glm_params.abs().sort_values(ascending=False).index[:20]
    glm_or_top = np.exp(glm_params.reindex(order))

    gain_dict = xgb_model.get_booster().get_score(importance_type="gain")
    # If the booster only has positional names (f0, f1, ...) map back to
    # real column names via all_columns' order.
    if gain_dict and all(k.startswith("f") and k[1:].isdigit() for k in gain_dict.keys()):
        idx_map   = {f"f{i}": col for i, col in enumerate(all_columns)}
        gain_dict = {idx_map.get(k, k): v for k, v in gain_dict.items()}
    gain_top = pd.Series(gain_dict).sort_values(ascending=False)[:20]

    fig, axes = plt.subplots(1, 2, figsize=(14, 8))

    axes[0].barh(glm_or_top.index[::-1], glm_or_top.values[::-1], color="#2b8cbe")
    axes[0].axvline(1.0, color="grey", linestyle="--", linewidth=1)
    axes[0].set_xlabel("Odds ratio")
    axes[0].set_title("GLM: top 20 features by |log-odds|")

    axes[1].barh(gain_top.index[::-1], gain_top.values[::-1], color="#de2d26")
    axes[1].set_xlabel("Gain importance")
    axes[1].set_title("XGBoost: top 20 features by gain")

    fig.tight_layout()

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "feature_importance_dual.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    print(f"[PLOT] Saved dual feature-importance plot -> {out_path}")
    return str(out_path)


def run_model_training_pipeline(selection_result: dict, cv_folds: int = 10, output_dir: str = "outputs") -> dict:
    """Calls compute_scale_pos_weight -> train_glm -> tune_xgboost ->
    cross_validated_metrics (both models) -> save_model_bundle ->
    plot_feature_importance_dual, in sequence.

    GLM trains on the ElasticNet-selected subset (interpretability).
    XGBoost trains on the FULL preprocessed feature set (predictive
    capacity) -- this is what all_train_columns in the saved bundle refers
    to, and what predict_cohort_dual()/imputer.transform() at inference
    time expect."""

    imputer            = selection_result["imputer"]
    scaler              = selection_result["scaler"]
    selected_features   = selection_result["selected_features"]
    X_train_scaled       = selection_result["X_train_scaled"]     # full columns -> XGBoost
    X_train_selected    = selection_result["X_train_selected"]    # reduced -> GLM
    y_train_raw          = selection_result["y_train"]

    y_binary = _to_binary(y_train_raw)
    y_binary.index = X_train_scaled.index

    scale_pos_weight = compute_scale_pos_weight(y_train_raw)

    glm = train_glm(X_train_selected, y_binary, selected_features, output_dir=output_dir)
    xgb_model = tune_xgboost(X_train_scaled, y_binary, scale_pos_weight, cv_folds=cv_folds)

    glm_cv = cross_validated_metrics(
        glm, X_train_selected, y_binary, cv_folds=cv_folds, model_name="glm", output_dir=output_dir
    )
    xgb_cv = cross_validated_metrics(
        xgb_model, X_train_scaled, y_binary, cv_folds=cv_folds, model_name="xgboost", output_dir=output_dir
    )

    all_train_columns = list(X_train_scaled.columns)
    training_metadata = {
        "n_train_samples":        int(X_train_scaled.shape[0]),
        "n_selected_features":    len(selected_features),
        "n_all_features":         len(all_train_columns),
        "feature_selection_seed": FEATURE_SELECTION_SEED,
        "xgboost_tuning_seed":    XGBOOST_TUNING_SEED,
        "scale_pos_weight":       scale_pos_weight,
        "glm_cv_auc_mean":        float(glm_cv["auc"].mean()),
        "glm_cv_auc_std":         float(glm_cv["auc"].std()),
        "xgb_cv_auc_mean":        float(xgb_cv["auc"].mean()),
        "xgb_cv_auc_std":         float(xgb_cv["auc"].std()),
    }

    bundle_path = save_model_bundle(
        glm=glm, xgb_model=xgb_model, imputer=imputer, scaler=scaler,
        selected_features=selected_features, all_train_columns=all_train_columns,
        training_metadata=training_metadata, output_dir="outputs/python",
    )

    plot_feature_importance_dual(
        glm, xgb_model, selected_features, all_train_columns, output_dir="outputs/python"
    )

    return {
        "glm":             glm,
        "xgb_model":       xgb_model,
        "glm_cv_results":  glm_cv,
        "xgb_cv_results":  xgb_cv,
        "bundle_path":     bundle_path,
    }