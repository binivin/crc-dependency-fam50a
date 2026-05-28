from pathlib import Path
import warnings
import numpy as np
import pandas as pd

from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


warnings.filterwarnings("ignore", category=ConvergenceWarning)

INTERIM_DIR = Path("data/interim")
TABLE_DIR = Path("artifacts/tables")
TABLE_DIR.mkdir(parents=True, exist_ok=True)

CANDIDATE_PATH = TABLE_DIR / "candidate_dependency_targets.csv"
TOP_N = 50


print("Reading data...")

cohort = pd.read_csv(INTERIM_DIR / "depmap_crc_cohort.csv")
gene_effect = pd.read_csv(INTERIM_DIR / "depmap_crc_gene_effect.csv")
expr = pd.read_csv(INTERIM_DIR / "depmap_crc_expression.csv")
cnv = pd.read_csv(INTERIM_DIR / "depmap_crc_cnv.csv")
mut = pd.read_csv(INTERIM_DIR / "depmap_crc_mutation_matrix.csv")
candidates = pd.read_csv(CANDIDATE_PATH)

for df in [cohort, gene_effect, expr, cnv, mut]:
    df["ModelID"] = df["ModelID"].astype(str)

# Top candidate targets
candidates = candidates.head(TOP_N).copy()
target_cols = candidates["target_column"].tolist()

print(f"Testing top {len(target_cols)} candidate targets")

# Prefix feature columns
expr_feat = expr.rename(columns={c: f"expr_{c}" for c in expr.columns if c != "ModelID"})
cnv_feat = cnv.rename(columns={c: f"cnv_{c}" for c in cnv.columns if c != "ModelID"})
mut_feat = mut.copy()

feature_data = (
    cohort[["ModelID", "PatientID", "CellLineName"]]
    .merge(expr_feat, on="ModelID", how="inner")
    .merge(cnv_feat, on="ModelID", how="inner")
    .merge(mut_feat, on="ModelID", how="inner")
)

feature_cols = [
    c for c in feature_data.columns
    if c not in ["ModelID", "PatientID", "CellLineName"]
]

X_all = feature_data[feature_cols].replace([np.inf, -np.inf], np.nan).fillna(0)
groups_all = feature_data["PatientID"].fillna(feature_data["ModelID"]).astype(str)

print(f"Samples: {X_all.shape[0]}")
print(f"Raw features: {X_all.shape[1]}")

cv = GroupKFold(n_splits=5)

# n=59라서 너무 많은 feature를 쓰지 않음
k = min(200, X_all.shape[1])

model = Pipeline([
    ("select", SelectKBest(score_func=f_regression, k=k)),
    ("scale", StandardScaler()),
    ("elastic", ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=50000, random_state=42)),
])

results = []

for i, target_col in enumerate(target_cols, start=1):
    gene_symbol = candidates.loc[candidates["target_column"] == target_col, "gene_symbol"].iloc[0]

    y_df = gene_effect[["ModelID", target_col]].rename(columns={target_col: "y"})

    data = (
        feature_data[["ModelID", "PatientID", "CellLineName"]]
        .merge(y_df, on="ModelID", how="inner")
    )

    y = data["y"].astype(float)
    valid = ~y.isna()

    X = X_all.loc[valid].copy()
    y = y.loc[valid].copy()
    groups = groups_all.loc[valid].copy()

    pred = cross_val_predict(model, X, y, cv=cv, groups=groups)

    mae = mean_absolute_error(y, pred)
    rmse = mean_squared_error(y, pred) ** 0.5
    r2 = r2_score(y, pred)

    baseline_pred = np.repeat(y.mean(), len(y))
    baseline_rmse = mean_squared_error(y, baseline_pred) ** 0.5

    results.append({
        "rank": i,
        "gene_symbol": gene_symbol,
        "target_column": target_col,
        "n_samples": len(y),
        "y_mean": y.mean(),
        "y_std": y.std(),
        "y_min": y.min(),
        "y_max": y.max(),
        "MAE": mae,
        "RMSE": rmse,
        "baseline_RMSE": baseline_rmse,
        "R2": r2,
    })

    print(f"[{i:02d}/{TOP_N}] {gene_symbol:10s} R2={r2: .4f}, RMSE={rmse:.4f}")

results_df = pd.DataFrame(results).sort_values("R2", ascending=False)

out_path = TABLE_DIR / "candidate_targets_elastic_net_results.csv"
results_df.to_csv(out_path, index=False, encoding="utf-8-sig")

print("\n=== Top results by R2 ===")
print(
    results_df[
        ["gene_symbol", "target_column", "y_std", "y_min", "y_max", "RMSE", "baseline_RMSE", "R2"]
    ].head(20).to_string(index=False)
)

print(f"\nSaved results to: {out_path}")