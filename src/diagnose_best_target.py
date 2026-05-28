from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.base import clone
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


warnings.filterwarnings("ignore", category=ConvergenceWarning)

INTERIM_DIR = Path("data/interim")
TABLE_DIR = Path("artifacts/tables")
FIG_DIR = Path("artifacts/figures")
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

RESULT_PATH = TABLE_DIR / "candidate_targets_elastic_net_results.csv"

print("Reading previous candidate target results...")
results = pd.read_csv(RESULT_PATH)

best = results.sort_values("R2", ascending=False).iloc[0]
TARGET_GENE = best["gene_symbol"]
TARGET_COL = best["target_column"]

print(f"Best target gene: {TARGET_GENE}")
print(f"Target column: {TARGET_COL}")
print(f"Previous R2: {best['R2']:.4f}")

print("\nReading CRC matrices...")

cohort = pd.read_csv(INTERIM_DIR / "depmap_crc_cohort.csv")
gene_effect = pd.read_csv(INTERIM_DIR / "depmap_crc_gene_effect.csv")
expr = pd.read_csv(INTERIM_DIR / "depmap_crc_expression.csv")
cnv = pd.read_csv(INTERIM_DIR / "depmap_crc_cnv.csv")
mut = pd.read_csv(INTERIM_DIR / "depmap_crc_mutation_matrix.csv")

for df in [cohort, gene_effect, expr, cnv, mut]:
    df["ModelID"] = df["ModelID"].astype(str)

# Prefix feature columns
expr_feat = expr.rename(columns={c: f"expr_{c}" for c in expr.columns if c != "ModelID"})
cnv_feat = cnv.rename(columns={c: f"cnv_{c}" for c in cnv.columns if c != "ModelID"})
mut_feat = mut.copy()

y_df = gene_effect[["ModelID", TARGET_COL]].rename(columns={TARGET_COL: "y"})

data = (
    cohort[["ModelID", "PatientID", "CellLineName"]]
    .merge(y_df, on="ModelID", how="inner")
    .merge(expr_feat, on="ModelID", how="inner")
    .merge(cnv_feat, on="ModelID", how="inner")
    .merge(mut_feat, on="ModelID", how="inner")
)

data = data.dropna(subset=["y"]).copy()

feature_cols = [
    c for c in data.columns
    if c not in ["ModelID", "PatientID", "CellLineName", "y"]
]

X = data[feature_cols].replace([np.inf, -np.inf], np.nan).fillna(0)
y = data["y"].astype(float)
groups = data["PatientID"].fillna(data["ModelID"]).astype(str)

print(f"Samples: {X.shape[0]}")
print(f"Raw features: {X.shape[1]}")
print(f"Target y range: {y.min():.3f} to {y.max():.3f}")
print(f"Target y std: {y.std():.3f}")

k = min(200, X.shape[1])

base_model = Pipeline([
    ("var", VarianceThreshold(threshold=0.0)),
    ("select", SelectKBest(score_func=f_regression, k=k)),
    ("scale", StandardScaler()),
    ("elastic", ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=50000, random_state=42)),
])

cv = GroupKFold(n_splits=5)

pred = np.zeros(len(y))
fold_records = []
feature_records = []

feature_array = np.array(feature_cols)

print("\nRunning fold-by-fold GroupKFold diagnostics...")

for fold, (train_idx, test_idx) in enumerate(cv.split(X, y, groups), start=1):
    model = clone(base_model)

    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

    model.fit(X_train, y_train)
    fold_pred = model.predict(X_test)
    pred[test_idx] = fold_pred

    fold_mae = mean_absolute_error(y_test, fold_pred)
    fold_rmse = mean_squared_error(y_test, fold_pred) ** 0.5

    # Fold R2 can be unstable with small test sets, so treat it cautiously.
    fold_r2 = r2_score(y_test, fold_pred)

    fold_records.append({
        "fold": fold,
        "n_train": len(train_idx),
        "n_test": len(test_idx),
        "MAE": fold_mae,
        "RMSE": fold_rmse,
        "R2": fold_r2,
    })

    var_mask = model.named_steps["var"].get_support()
    after_var_features = feature_array[var_mask]

    select_mask = model.named_steps["select"].get_support()
    selected_features = after_var_features[select_mask]

    coefs = model.named_steps["elastic"].coef_

    for feat, coef in zip(selected_features, coefs):
        if abs(coef) > 1e-8:
            feature_records.append({
                "fold": fold,
                "feature": feat,
                "coef": coef,
                "abs_coef": abs(coef),
            })

    print(f"Fold {fold}: RMSE={fold_rmse:.4f}, R2={fold_r2:.4f}")

mae = mean_absolute_error(y, pred)
rmse = mean_squared_error(y, pred) ** 0.5
r2 = r2_score(y, pred)

baseline_pred = np.repeat(y.mean(), len(y))
baseline_rmse = mean_squared_error(y, baseline_pred) ** 0.5

corr = np.corrcoef(y, pred)[0, 1]

print("\n=== Overall diagnostics ===")
print(f"Target: {TARGET_GENE}")
print(f"MAE: {mae:.4f}")
print(f"RMSE: {rmse:.4f}")
print(f"Baseline RMSE: {baseline_rmse:.4f}")
print(f"R2: {r2:.4f}")
print(f"Correlation: {corr:.4f}")

# Save prediction table
pred_df = data[["ModelID", "PatientID", "CellLineName", "y"]].copy()
pred_df["pred"] = pred
pred_df["residual"] = pred_df["y"] - pred_df["pred"]

pred_path = TABLE_DIR / f"diagnostic_{TARGET_GENE}_predictions.csv"
pred_df.to_csv(pred_path, index=False, encoding="utf-8-sig")

# Save fold metrics
fold_df = pd.DataFrame(fold_records)
fold_path = TABLE_DIR / f"diagnostic_{TARGET_GENE}_fold_metrics.csv"
fold_df.to_csv(fold_path, index=False, encoding="utf-8-sig")

# Save feature stability
feature_df = pd.DataFrame(feature_records)

if len(feature_df) > 0:
    stability = (
        feature_df
        .groupby("feature")
        .agg(
            n_folds=("fold", "nunique"),
            mean_abs_coef=("abs_coef", "mean"),
            mean_coef=("coef", "mean")
        )
        .reset_index()
        .sort_values(["n_folds", "mean_abs_coef"], ascending=False)
    )

    stability_path = TABLE_DIR / f"diagnostic_{TARGET_GENE}_feature_stability.csv"
    stability.to_csv(stability_path, index=False, encoding="utf-8-sig")

    print("\nTop stable features:")
    print(stability.head(20).to_string(index=False))
else:
    print("\nNo nonzero features selected.")

# Plot actual vs predicted
plt.figure(figsize=(6, 6))
plt.scatter(y, pred)

min_val = min(y.min(), pred.min())
max_val = max(y.max(), pred.max())
plt.plot([min_val, max_val], [min_val, max_val], linestyle="--")

plt.xlabel("Actual CRISPR gene effect")
plt.ylabel("Predicted CRISPR gene effect")
plt.title(f"{TARGET_GENE} dependency prediction\nR2={r2:.3f}, corr={corr:.3f}")
plt.tight_layout()

fig_path = FIG_DIR / f"diagnostic_{TARGET_GENE}_actual_vs_predicted.png"
plt.savefig(fig_path, dpi=300)

print(f"\nSaved predictions to: {pred_path}")
print(f"Saved fold metrics to: {fold_path}")
print(f"Saved figure to: {fig_path}")