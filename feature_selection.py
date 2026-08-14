
from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.feature_selection import mutual_info_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler


def fit_preprocessing(X_train: pd.DataFrame):
    """Fit SimpleImputer(strategy="median") and StandardScaler() on
    X_train ONLY -- never on test -- so no test-set statistics leak into
    preprocessing. Returns the fitted (imputer, scaler) pair."""
    imputer = SimpleImputer(strategy="median")
    X_train_imputed = imputer.fit_transform(X_train)

    scaler = StandardScaler()
    scaler.fit(X_train_imputed)

    return imputer, scaler


def apply_preprocessing(X: pd.DataFrame, imputer, scaler) -> pd.DataFrame:
    """Transform-only -- never fit. Applies preprocessing fitted elsewhere
    (i.e. on train) to any matrix, including test, without refitting."""
    X_imputed = imputer.transform(X)
    X_scaled  = scaler.transform(X_imputed)
    return pd.DataFrame(X_scaled, index=X.index, columns=X.columns)


def run_elasticnet_selection(X_train_scaled: pd.DataFrame, y_train: pd.Series,
                              cv_folds: int = 10, random_state: int = 42,
                              output_dir: str = "outputs") -> dict:
    """ElasticNetCV used purely for feature selection (exact-zero
    coefficients), not as a predictive model -- the GLM and XGBoost in
    model_training.py are the actual predictors."""

    y_aligned = pd.Series(y_train).reindex(X_train_scaled.index) if hasattr(y_train, "reindex") \
        else pd.Series(y_train, index=X_train_scaled.index)
    n_unmatched = int(y_aligned.isna().sum())
    if n_unmatched:
        raise ValueError(
            f"[ELASTICNET] {n_unmatched} sample(s) in X_train_scaled have no matching label after "
            f"alignment; refusing to train with an incomplete y. Check that y_train is keyed on the "
            f"same sample_id index as X_train_scaled."
        )
    y_binary = y_aligned.astype(str).str.strip().eq("pCR").astype(int)

    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=random_state)

    model = ElasticNetCV(
        l1_ratio=[0.1, 0.3, 0.5, 0.7, 0.9],
        cv=cv,
        max_iter=10000,
        n_jobs=-1,
        random_state=random_state,
    )
    # Threading backend: avoids spawning loky worker processes (each of which
    # re-imports the full stack under macOS `spawn`) and the CPython
    # resource_tracker shutdown race that raised `ChildProcessError` at exit.
    with joblib.parallel_config(backend="threading", n_jobs=-1):
        model.fit(X_train_scaled.values, y_binary.values)

    coef_table = pd.DataFrame({
        "feature": X_train_scaled.columns,
        "coefficient": model.coef_,
    })
    coef_table["selected"] = coef_table["coefficient"] != 0
    coef_table = coef_table.reindex(
        coef_table["coefficient"].abs().sort_values(ascending=False).index
    ).reset_index(drop=True)

    selected_features = coef_table.loc[coef_table["selected"], "feature"].tolist()
    if not selected_features:
        raise ValueError(
            "[ELASTICNET] ElasticNet zeroed every coefficient; no features selected. "
            "Relax the l1_ratio range or check X/y alignment."
        )

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "selected_features.json", "w") as f:
        json.dump(selected_features, f, indent=2)
    coef_table.to_csv(out_dir / "elasticnet_coefficients.csv", index=False)

    print(
        f"[ELASTICNET] Retained {len(selected_features)} / {X_train_scaled.shape[1]} features "
        f"(l1_ratio={model.l1_ratio_:.2f}, alpha={model.alpha_:.5f})."
    )

    return {
        "selected_features": selected_features,
        "coefficients": coef_table,
        "model": model,
    }


def run_feature_selection_pipeline(fused_train: pd.DataFrame, labels_train: pd.Series,
                                    cv_folds: int = 10, random_state: int = 42,
                                    output_dir: str = "outputs",
                                    drop_gse_id: bool = False,
                                    feature_reduction_k: int | None = None) -> dict:
    working = fused_train
    if drop_gse_id:
        gse_cols = [c for c in working.columns if c.startswith("gse_id_")]
        if gse_cols:
            working = working.drop(columns=gse_cols)
            print(f"[SELECTION] Dropped {len(gse_cols)} gse_id (study) columns for cohort-agnostic modeling.")

    if feature_reduction_k is not None:
        # Impute/scale once (on `working`) only to compute MI scores, then
        # keep the top-k columns and FIT the shipped preprocessing on that
        # reduced set so imputer/scaler stay consistent with the columns the
        # models (and inference) actually use.
        tmp_imputer, tmp_scaler = fit_preprocessing(working)
        X_tmp = apply_preprocessing(working, tmp_imputer, tmp_scaler)
        y_aligned = pd.Series(labels_train).reindex(working.index) if hasattr(labels_train, "reindex") \
            else pd.Series(labels_train, index=working.index)
        y_binary = y_aligned.astype(str).str.strip().eq("pCR").astype(int)
        mi = mutual_info_classif(
            X_tmp.values, y_binary.values, random_state=random_state, n_neighbors=3
        )
        keep = sorted(np.argsort(mi)[-feature_reduction_k:])
        keep_cols = [X_tmp.columns[i] for i in keep]
        working = working[keep_cols]
        print(
            f"[SELECTION] Mutual-information reduction: kept top {len(keep_cols)} / "
            f"{X_tmp.shape[1]} features."
        )

    imputer, scaler = fit_preprocessing(working)
    X_train_scaled = apply_preprocessing(working, imputer, scaler)

    selection = run_elasticnet_selection(
        X_train_scaled, labels_train, cv_folds=cv_folds,
        random_state=random_state, output_dir=output_dir
    )

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    preproc_path = out_dir / "preprocessing.joblib"
    joblib.dump({"imputer": imputer, "scaler": scaler}, preproc_path)
    print(f"[PREPROCESS] Saved {preproc_path}")

    X_train_selected = X_train_scaled[selection["selected_features"]]
    y_train_aligned = pd.Series(labels_train).reindex(fused_train.index) if hasattr(labels_train, "reindex") \
        else pd.Series(labels_train, index=fused_train.index)

    return {
        "imputer": imputer,
        "scaler": scaler,
        "selected_features": selection["selected_features"],
        "X_train_scaled": X_train_scaled,
        "X_train_selected": X_train_selected,
        "y_train": y_train_aligned,
    }