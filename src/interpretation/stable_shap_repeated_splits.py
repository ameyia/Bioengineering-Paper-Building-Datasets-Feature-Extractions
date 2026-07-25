import argparse
import os
import warnings
from collections import Counter

import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    matthews_corrcoef,
    accuracy_score,
    f1_score,
    precision_score,
    recall_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier

warnings.filterwarnings("ignore")

try:
    import shap
except ImportError:
    raise ImportError("Please install SHAP first: pip install shap")

try:
    from xgboost import XGBClassifier
    HAS_XGB = True
except ImportError:
    HAS_XGB = False


def get_models(random_state):
    models = {
        "random_forest": RandomForestClassifier(
            n_estimators=500,
            max_depth=None,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=random_state,
            n_jobs=-1,
        ),
        "logistic_regression": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        max_iter=5000,
                        class_weight="balanced",
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "svm_rbf": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    SVC(
                        kernel="rbf",
                        C=10,
                        gamma="scale",
                        probability=True,
                        class_weight="balanced",
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "knn": Pipeline(
            [
                ("scaler", StandardScaler()),
                ("model", KNeighborsClassifier(n_neighbors=5)),
            ]
        ),
    }

    if HAS_XGB:
        models["xgboost"] = XGBClassifier(
            n_estimators=300,
            max_depth=3,
            learning_rate=0.03,
            subsample=0.8,
            colsample_bytree=0.8,
            eval_metric="logloss",
            random_state=random_state,
            n_jobs=-1,
        )

    return models


def evaluate_model(model, X, y):
    pred = model.predict(X)

    return {
        "mcc": matthews_corrcoef(y, pred),
        "accuracy": accuracy_score(y, pred),
        "f1": f1_score(y, pred, zero_division=0),
        "precision": precision_score(y, pred, zero_division=0),
        "recall": recall_score(y, pred, zero_division=0),
    }


def get_positive_class_shap_values(shap_values):
    """
    Handles different SHAP output formats:
    - list[class0, class1]
    - array n_samples x n_features
    - array n_samples x n_features x classes
    """
    if isinstance(shap_values, list):
        return shap_values[1]

    shap_values = np.array(shap_values)

    if shap_values.ndim == 3:
        return shap_values[:, :, 1]

    return shap_values


def compute_shap_importance(model, X_train, X_test, feature_names, max_shap_samples=100):
    """
    Uses TreeExplainer for tree models.
    Uses KernelExplainer for non-tree models, but limits samples because it is slower.
    """

    X_test_sample = X_test.copy()

    if len(X_test_sample) > max_shap_samples:
        X_test_sample = X_test_sample.sample(
            n=max_shap_samples,
            random_state=1,
        )

    # If pipeline, use the final estimator type to decide.
    final_model = model
    transformed_train = X_train
    transformed_test = X_test_sample

    if isinstance(model, Pipeline):
        final_model = model.named_steps["model"]

        if "scaler" in model.named_steps:
            transformed_train = model.named_steps["scaler"].transform(X_train)
            transformed_test = model.named_steps["scaler"].transform(X_test_sample)

    tree_model_names = [
        "RandomForestClassifier",
        "XGBClassifier",
        "DecisionTreeClassifier",
        "ExtraTreesClassifier",
        "GradientBoostingClassifier",
    ]

    if final_model.__class__.__name__ in tree_model_names:
        explainer = shap.TreeExplainer(final_model)
        shap_values = explainer.shap_values(transformed_test)
        shap_class1 = get_positive_class_shap_values(shap_values)

    elif final_model.__class__.__name__ == "LogisticRegression":
        explainer = shap.LinearExplainer(final_model, transformed_train)
        shap_values = explainer.shap_values(transformed_test)
        shap_class1 = get_positive_class_shap_values(shap_values)

    else:
        # Slower fallback for SVM/KNN
        background_for_shap = X_train.sample(
        n=min(50, len(X_train)),
        random_state=1,
        )

        test_for_shap = X_test_sample.copy()

        # Wrap predict_proba so SHAP does not directly modify sklearn attributes
        def predict_proba_wrapper(data):
            data_df = pd.DataFrame(data, columns=feature_names)
            return model.predict_proba(data_df)

        explainer = shap.KernelExplainer(
            predict_proba_wrapper,
            background_for_shap,
        )

        shap_values = explainer.shap_values(
            test_for_shap,
            nsamples=100,
        )

        shap_class1 = get_positive_class_shap_values(shap_values)

    mean_abs_shap = np.abs(shap_class1).mean(axis=0)

    importance = pd.DataFrame(
        {
            "feature": feature_names,
            "mean_abs_shap": mean_abs_shap,
        }
    ).sort_values("mean_abs_shap", ascending=False)

    return importance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True, help="Feature table TSV file")
    parser.add_argument("--outdir", default="results/interpretability/stable_shap_results")
    parser.add_argument("--n_splits", type=int, default=10)
    parser.add_argument("--top_n", type=int, default=20)
    parser.add_argument("--max_shap_samples", type=int, default=100)

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.table, sep="\t")

    if "label" not in df.columns:
        raise ValueError("The table must contain a 'label' column.")

    y = df["label"].astype(int)

    drop_cols = ["label"]
    if "sequence" in df.columns:
        drop_cols.append("sequence")

    X = df.drop(columns=drop_cols)

    # Make sure all features are numeric.
    X = X.apply(pd.to_numeric, errors="coerce")
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median(numeric_only=True))

    feature_names = X.columns.tolist()

    split_metrics = []
    all_top_features = []
    all_shap_rows = []

    seeds = list(range(args.n_splits))

    for seed in seeds:
        print(f"\n===== Split seed {seed} =====")

        # First split: 60% train, 40% temporary
        X_train, X_temp, y_train, y_temp = train_test_split(
            X,
            y,
            test_size=0.40,
            stratify=y,
            random_state=seed,
        )

        # Second split: 20% validation, 20% test
        X_val, X_test, y_val, y_test = train_test_split(
            X_temp,
            y_temp,
            test_size=0.50,
            stratify=y_temp,
            random_state=seed,
        )

        models = get_models(seed)

        val_results = []

        for model_name, model in models.items():
            model.fit(X_train, y_train)

            val_metrics = evaluate_model(model, X_val, y_val)
            val_results.append(
                {
                    "seed": seed,
                    "model": model_name,
                    **{f"val_{k}": v for k, v in val_metrics.items()},
                }
            )

        val_results_df = pd.DataFrame(val_results)
        best_row = val_results_df.sort_values(
            "val_mcc",
            ascending=False,
        ).iloc[0]

        best_model_name = best_row["model"]
        best_model = models[best_model_name]

        test_metrics = evaluate_model(best_model, X_test, y_test)

        print(f"Best model: {best_model_name}")
        print(f"Validation MCC: {best_row['val_mcc']:.3f}")
        print(f"Test MCC: {test_metrics['mcc']:.3f}")

        split_metrics.append(
            {
                "seed": seed,
                "best_model": best_model_name,
                **best_row.to_dict(),
                **{f"test_{k}": v for k, v in test_metrics.items()},
            }
        )

        shap_importance = compute_shap_importance(
            best_model,
            X_train,
            X_test,
            feature_names,
            max_shap_samples=args.max_shap_samples,
        )

        shap_importance["seed"] = seed
        shap_importance["best_model"] = best_model_name
        shap_importance["rank"] = np.arange(1, len(shap_importance) + 1)

        top_features = shap_importance.head(args.top_n).copy()
        all_top_features.extend(top_features["feature"].tolist())
        all_shap_rows.append(shap_importance)

        top_features.to_csv(
            os.path.join(args.outdir, f"top_shap_features_seed_{seed}.csv"),
            index=False,
        )

    metrics_df = pd.DataFrame(split_metrics)
    metrics_df.to_csv(
        os.path.join(args.outdir, "split_metrics.csv"),
        index=False,
    )

    all_shap_df = pd.concat(all_shap_rows, ignore_index=True)
    all_shap_df.to_csv(
        os.path.join(args.outdir, "all_shap_importance_by_split.csv"),
        index=False,
    )

    # Stable feature summary
    top_counter = Counter(all_top_features)

    stable_rows = []

    for feature, count in top_counter.items():
        feature_rows = all_shap_df[all_shap_df["feature"] == feature]

        stable_rows.append(
            {
                "feature": feature,
                f"times_in_top_{args.top_n}": count,
                "fraction_of_splits": count / args.n_splits,
                "mean_rank": feature_rows["rank"].mean(),
                "median_rank": feature_rows["rank"].median(),
                "mean_abs_shap_across_splits": feature_rows["mean_abs_shap"].mean(),
                "std_abs_shap_across_splits": feature_rows["mean_abs_shap"].std(),
            }
        )

    stable_df = pd.DataFrame(stable_rows).sort_values(
        [f"times_in_top_{args.top_n}", "mean_abs_shap_across_splits"],
        ascending=[False, False],
    )

    stable_df.to_csv(
        os.path.join(args.outdir, "stable_shap_features.csv"),
        index=False,
    )

    print("\nDone.")
    print(f"Results saved to: {args.outdir}")
    print("\nMost stable SHAP features:")
    print(stable_df.head(20))


if __name__ == "__main__":
    main()