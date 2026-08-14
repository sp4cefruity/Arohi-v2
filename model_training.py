from __future__ import annotations

from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import optuna
from sklearn.base import clone as sk_clone
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from xgboost import XGBClassifier

FEATURE_SELECTION_SEED = 42   # ElasticNet selection / pre-processing
XGBOOST_TUNING_SEED     = 7    # XGBoost Optuna search
REPORTING_SEED          = 123  # CV folds for the reported metrics -- must NOT
                               # overlap the tuning folds, or the reported AUC
                               # is optimistically biased (params picked on the
                               # very folds being scored).


def compute_scale_pos_weight(y_train) -> float:
    """n_RD / n_pCR, used as the centre of XGBoost's scale_pos_weight search
    range. Accepts either string labels ('pCR'/'RD') or a binary encoding
    (1=pCR, 0=RD)."""
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


def _to_binary(y) -> pd.Series:
    """pCR -> 1, RD -> 0. Passes through an already-binary series unchanged."""
    s = pd.Series(y)
    if s.dtype == object or str(s.dtype).startswith("string") or str(s.dtype) == "category":
        return (s.astype(str).str.strip() == "pCR").astype(int)
    return s.astype(int)


def train_glm(X_train: pd.DataFrame, y_binary, selected_features: list[str],
              output_dir: str = "outputs"):
    """Regularized logistic regression on the ElasticNet-selected features.

    sklearn's L2 LogisticRegression replaces the statsmodels Logit used
    before, which repeatedly failed to converge within maxiter=500 on this
    near-separated feature set. C is tuned with a quick 5-fold CV. Odds
    ratios stay interpretable (exp(coef))."""
    X = X_train[selected_features]

    Cs = np.logspace(-3, 2, 12)
    best_c, best_auc = Cs[0], -1.0
    c_cv = StratifiedKFold(5, shuffle=True, random_state=FEATURE_SELECTION_SEED)
    for c in Cs:
        probe = LogisticRegression(
            penalty="l2", C=c, class_weight="balanced",
            solver="lbfgs", max_iter=2000, random_state=FEATURE_SELECTION_SEED,
        )
        scores = cross_val_score(probe, X, y_binary, cv=c_cv, scoring="roc_auc", n_jobs=1)
        if scores.mean() > best_auc:
            best_auc, best_c = float(scores.mean()), c

    model = LogisticRegression(
        penalty="l2", C=best_c, class_weight="balanced",
        solver="lbfgs", max_iter=2000, random_state=FEATURE_SELECTION_SEED,
    )
    model.fit(X, y_binary)

    coef_table = pd.DataFrame({
        "feature":     selected_features,
        "coefficient": model.coef_[0],
    })
    coef_table["odds_ratio"] = np.exp(coef_table["coefficient"])

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "glm_coefficients.csv"
    coef_table.to_csv(out_path, index=False)

    print(
        f"[GLM] Fit regularized logistic regression on {len(selected_features)} selected feature(s) "
        f"(C={best_c:.4f}, 5-fold CV AUC={best_auc:.4f}) -> {out_path}"
    )

    return model


def tune_xgboost(X_train: pd.DataFrame, y_binary, scale_pos_weight: float,
                 cv_folds: int = 5, n_trials: int = 60,
                 seed: int = XGBOOST_TUNING_SEED):
    """Optuna search over the XGBoost hyperparameters with per-fold early
    stopping (n_estimators self-caps via best_iteration). Refits the best
    config on all training data afterwards. Returns the fitted model and the
    study object (for best params/score)."""
    y = np.asarray(pd.Series(y_binary))
    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=seed)

    def objective(trial):
        params = {
            "objective": "binary:logistic",
            "eval_metric": "auc",
            "tree_method": "hist",
            "n_estimators": 2000,
            "early_stopping_rounds": 30,
            "n_jobs": 1,
            "random_state": seed,
            "max_depth":        trial.suggest_int("max_depth", 2, 6),
            "learning_rate":    trial.suggest_float("learning_rate", 0.005, 0.15, log=True),
            "subsample":        trial.suggest_float("subsample", 0.5, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.3, 0.9),
            "colsample_bylevel": trial.suggest_float("colsample_bylevel", 0.5, 1.0),
            "min_child_weight": trial.suggest_int("min_child_weight", 1, 10),
            "gamma":            trial.suggest_float("gamma", 0.0, 1.0),
            "reg_alpha":        trial.suggest_float("reg_alpha", 1e-3, 10.0, log=True),
            "reg_lambda":       trial.suggest_float("reg_lambda", 1e-3, 10.0, log=True),
            "max_delta_step":   trial.suggest_int("max_delta_step", 0, 5),
            "scale_pos_weight": trial.suggest_float(
                "scale_pos_weight", 0.5 * scale_pos_weight, 1.5 * scale_pos_weight
            ),
        }
        aucs, best_iters = [], []
        for train_idx, val_idx in cv.split(X_train, y):
            X_tr, X_va = X_train.iloc[train_idx], X_train.iloc[val_idx]
            y_tr, y_va = y[train_idx], y[val_idx]
            model = XGBClassifier(**params)
            model.fit(
                X_tr, y_tr,
                eval_set=[(X_va, y_va)],
                verbose=False,
            )
            aucs.append(roc_auc_score(y_va, model.predict_proba(X_va)[:, 1]))
            if getattr(model, "best_iteration", None) is not None:
                best_iters.append(int(model.best_iteration) + 1)
        trial.set_user_attr("best_iter", int(np.mean(best_iters)) if best_iters else 2000)
        return float(np.mean(aucs))

    sampler = optuna.samplers.TPESampler(seed=seed)
    study = optuna.create_study(direction="maximize", sampler=sampler)
    study.optimize(objective, n_trials=n_trials, show_progress_bar=False)

    best_params = dict(study.best_params)
    best_iter = int(study.best_trial.user_attrs.get("best_iter", 2000))
    best_params["n_estimators"] = best_iter

    final = XGBClassifier(
        objective="binary:logistic",
        eval_metric="auc",
        tree_method="hist",
        n_jobs=1,
        random_state=seed,
        **best_params,
    )
    final.fit(X_train, y)

    print(f"[XGB TUNE] Best CV AUC: {study.best_value:.4f}")
    print(f"[XGB TUNE] Best params: {best_params}")

    return final, study


def cross_validated_metrics(model, X_train: pd.DataFrame, y_binary, cv_folds: int = 10,
                             model_name: str = "model", output_dir: str = "outputs",
                             seed: int = REPORTING_SEED) -> pd.DataFrame:
    """Out-of-fold AUC / AUPRC / F1 / accuracy. Refits a fresh clone of
    `model` inside each fold so metrics aren't leaked from the already-fit
    model passed in. Uses REPORTING_SEED folds -- deliberately different
    from the tuning seed (and from feature selection's 42)."""
    y = np.asarray(pd.Series(y_binary))
    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=seed)

    rows = []
    for fold_idx, (train_idx, val_idx) in enumerate(cv.split(X_train, y), start=1):
        X_fold_train, X_fold_val = X_train.iloc[train_idx], X_train.iloc[val_idx]
        y_fold_train, y_fold_val = y[train_idx], y[val_idx]

        fold_model = sk_clone(model)
        fold_model.fit(X_fold_train, y_fold_train)
        y_pred_proba = fold_model.predict_proba(X_fold_val)[:, 1]
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


def generate_oof_probabilities(model, X_train: pd.DataFrame, y_binary, cv_folds: int = 5,
                                seed: int = XGBOOST_TUNING_SEED) -> pd.Series:
    """Out-of-fold predicted probabilities from a fresh clone of `model`.
    X_train and y_binary must be index-aligned (both sliced by the same
    row subset if used inside a nested loop)."""
    y = pd.Series(np.asarray(pd.Series(y_binary)), index=X_train.index)
    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=seed)
    oof = pd.Series(np.nan, index=X_train.index, dtype=float)

    for train_idx, val_idx in cv.split(X_train, y.values):
        X_fold_train, X_fold_val = X_train.iloc[train_idx], X_train.iloc[val_idx]
        fold_model = sk_clone(model)
        fold_model.fit(X_fold_train, y.iloc[train_idx])
        oof.iloc[val_idx] = fold_model.predict_proba(X_fold_val)[:, 1]

    return oof


def train_stacking_meta(oof_probs: pd.DataFrame, y_binary) -> LogisticRegression:
    """Meta-learner: logistic regression on the OOF probabilities of the base
    models (GLM, XGBoost). Fitted strictly on out-of-fold predictions so its
    training signal is not leaked into the base-model scores."""
    y = np.asarray(pd.Series(y_binary)).astype(int)
    meta = LogisticRegression(penalty="l2", C=1.0, max_iter=1000)
    meta.fit(oof_probs, y)
    print(f"[STACK] Meta-learner (LR on {oof_probs.shape[1]} OOF prob column(s)) fitted.")
    return meta


def cross_validated_stack_metrics(glm, xgb_model, X_train_selected: pd.DataFrame,
                                   X_train_scaled: pd.DataFrame, y_binary,
                                   cv_folds: int = 10, inner_folds: int = 5,
                                   output_dir: str = "outputs") -> pd.DataFrame:
    """Honest nested-CV metrics for the stacking ensemble. For every outer
    fold: the base models are re-fit on the outer train, their OOF
    probabilities on the inner folds train the meta-learner, and the outer
    validation is scored by the full stack. No fold shared between the
    tuning folds (7), the inner OOF folds (7) and the reporting folds (123)
    leaks the meta-model's training signal into its score."""
    y = np.asarray(pd.Series(y_binary))
    outer = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=REPORTING_SEED)

    rows = []
    for fold_idx, (tr, va) in enumerate(outer.split(X_train_scaled, y), start=1):
        X_tr_sel, X_va_sel = X_train_selected.iloc[tr], X_train_selected.iloc[va]
        X_tr_all, X_va_all = X_train_scaled.iloc[tr], X_train_scaled.iloc[va]
        y_tr, y_va = y[tr], y[va]
        y_tr_series = pd.Series(y_tr, index=X_tr_all.index)

        glm_oof = generate_oof_probabilities(glm, X_tr_sel, y_tr_series, cv_folds=inner_folds)
        xgb_oof = generate_oof_probabilities(xgb_model, X_tr_all, y_tr_series, cv_folds=inner_folds)
        oof_inner = pd.DataFrame({"glm": glm_oof, "xgb": xgb_oof})
        meta = train_stacking_meta(oof_inner, y_tr_series)

        g = sk_clone(glm).fit(X_tr_sel, y_tr)
        x = sk_clone(xgb_model).fit(X_tr_all, y_tr)
        stack_proba = meta.predict_proba(pd.DataFrame({
            "glm": g.predict_proba(X_va_sel)[:, 1],
            "xgb": x.predict_proba(X_va_all)[:, 1],
        }))[:, 1]

        y_pred_label = (stack_proba >= 0.5).astype(int)
        rows.append({
            "fold":     fold_idx,
            "auc":      roc_auc_score(y_va, stack_proba),
            "auprc":    average_precision_score(y_va, stack_proba),
            "f1":       f1_score(y_va, y_pred_label, zero_division=0),
            "accuracy": accuracy_score(y_va, y_pred_label),
        })

    results = pd.DataFrame(rows)
    summary = results[["auc", "auprc", "f1", "accuracy"]].agg(["mean", "std"])
    print(
        f"[CV METRICS][stack] "
        f"AUC={summary.loc['mean','auc']:.3f}±{summary.loc['std','auc']:.3f}  "
        f"AUPRC={summary.loc['mean','auprc']:.3f}±{summary.loc['std','auprc']:.3f}  "
        f"F1={summary.loc['mean','f1']:.3f}±{summary.loc['std','f1']:.3f}  "
        f"Accuracy={summary.loc['mean','accuracy']:.3f}±{summary.loc['std','accuracy']:.3f}"
    )

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "cv_results_stack.csv"
    results.to_csv(out_path, index=False)
    print(f"[CV METRICS][stack] Saved -> {out_path}")

    return results


def save_model_bundle(glm, xgb_model, imputer, scaler, selected_features, all_train_columns,
                       stacking_meta, training_metadata: dict,
                       output_dir: str = "outputs/python") -> str:
    """Single joblib dump containing everything needed to reproduce
    predictions with no refitting."""
    bundle = {
        "glm":                glm,
        "xgb_model":          xgb_model,
        "imputer":            imputer,
        "scaler":             scaler,
        "selected_features":  selected_features,
        "all_train_columns":  all_train_columns,
        "stacking_meta":      stacking_meta,
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

    glm_coef = pd.Series(glm.coef_[0], index=selected_features)
    order    = glm_coef.abs().sort_values(ascending=False).index[:20]
    glm_or_top = np.exp(glm_coef.reindex(order))

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


def run_model_training_pipeline(selection_result: dict, cv_folds: int = 10,
                                output_dir: str = "outputs",
                                tune_cv_folds: int = 5, report_cv_folds: int = 10,
                                n_trials: int = 60,
                                drop_gse_id: bool = False,
                                feature_reduction_k: int | None = None) -> dict:
    """compute_scale_pos_weight -> train_glm -> tune_xgboost (Optuna) ->
    cross_validated_metrics (GLM + XGBoost, honest folds) -> OOF-probability
    stacking -> save_model_bundle -> plot_feature_importance_dual.

    GLM trains on the ElasticNet-selected subset (interpretability).
    XGBoost trains on the FULL preprocessed feature set (predictive
    capacity) -- this is what all_train_columns in the saved bundle refers
    to, and what predict_cohort_dual()/imputer.transform() at inference
    time expect."""

    imputer            = selection_result["imputer"]
    scaler              = selection_result["scaler"]
    selected_features   = selection_result["selected_features"]
    X_train_scaled       = selection_result["X_train_scaled"]
    X_train_selected    = selection_result["X_train_selected"]
    y_train_raw          = selection_result["y_train"]

    y_binary = _to_binary(y_train_raw)
    y_binary.index = X_train_scaled.index

    scale_pos_weight = compute_scale_pos_weight(y_train_raw)

    glm = train_glm(X_train_selected, y_binary, selected_features, output_dir=output_dir)
    xgb_model, study = tune_xgboost(
        X_train_scaled, y_binary, scale_pos_weight,
        cv_folds=tune_cv_folds, n_trials=n_trials,
    )

    glm_cv = cross_validated_metrics(
        glm, X_train_selected, y_binary, cv_folds=report_cv_folds,
        model_name="glm", output_dir=output_dir,
    )
    xgb_cv = cross_validated_metrics(
        xgb_model, X_train_scaled, y_binary, cv_folds=report_cv_folds,
        model_name="xgboost", output_dir=output_dir,
    )

    glm_oof = generate_oof_probabilities(glm, X_train_selected, y_binary)
    xgb_oof = generate_oof_probabilities(xgb_model, X_train_scaled, y_binary)
    oof_probs = pd.DataFrame({"glm": glm_oof, "xgb": xgb_oof}, index=X_train_scaled.index)
    stacking_meta = train_stacking_meta(oof_probs, y_binary)

    stack_cv = cross_validated_stack_metrics(
        glm, xgb_model, X_train_selected, X_train_scaled, y_binary,
        cv_folds=report_cv_folds, inner_folds=tune_cv_folds, output_dir=output_dir,
    )

    all_train_columns = list(X_train_scaled.columns)
    training_metadata = {
        "n_train_samples":        int(X_train_scaled.shape[0]),
        "n_selected_features":    len(selected_features),
        "n_all_features":         len(all_train_columns),
        "feature_selection_seed": FEATURE_SELECTION_SEED,
        "xgboost_tuning_seed":    XGBOOST_TUNING_SEED,
        "reporting_seed":         REPORTING_SEED,
        "scale_pos_weight":       scale_pos_weight,
        "drop_gse_id":            drop_gse_id,
        "feature_reduction_k":    feature_reduction_k,
        "xgb_best_params":        dict(study.best_params),
        "xgb_best_cv_auc":        float(study.best_value),
        "glm_cv_auc_mean":        float(glm_cv["auc"].mean()),
        "glm_cv_auc_std":         float(glm_cv["auc"].std()),
        "xgb_cv_auc_mean":        float(xgb_cv["auc"].mean()),
        "xgb_cv_auc_std":         float(xgb_cv["auc"].std()),
        "stack_cv_auc_mean":      float(stack_cv["auc"].mean()),
        "stack_cv_auc_std":       float(stack_cv["auc"].std()),
    }

    bundle_path = save_model_bundle(
        glm=glm, xgb_model=xgb_model, imputer=imputer, scaler=scaler,
        selected_features=selected_features, all_train_columns=all_train_columns,
        stacking_meta=stacking_meta, training_metadata=training_metadata,
        output_dir="outputs/python",
    )

    plot_feature_importance_dual(
        glm, xgb_model, selected_features, all_train_columns, output_dir="outputs/python"
    )

    return {
        "glm":             glm,
        "xgb_model":       xgb_model,
        "stacking_meta":   stacking_meta,
        "glm_cv_results":  glm_cv,
        "xgb_cv_results":  xgb_cv,
        "stack_cv_results": stack_cv,
        "oof_probs":       oof_probs,
        "training_metadata": training_metadata,
        "bundle_path":     bundle_path,
    }
