from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


INTERIM_DIR = Path("data/interim")
OUT_DIR = Path("artifacts/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)

TARGET_GENES = ["KRAS", "BRAF", "PIK3CA", "APC", "TP53", "SMAD4", "FBXW7"]


def find_gene_column(columns, gene):
    hits = []
    for col in columns:
        col = str(col)
        if col == gene or col.startswith(gene + " ") or col.startswith(gene + " ("):
            hits.append(col)
    return hits[0] if hits else None


print("Reading CRC matrices...")

cohort = pd.read_csv(INTERIM_DIR / "depmap_crc_cohort.csv")
gene_effect = pd.read_csv(INTERIM_DIR / "depmap_crc_gene_effect.csv")
expr = pd.read_csv(INTERIM_DIR / "depmap_crc_expression.csv")
cnv = pd.read_csv(INTERIM_DIR / "depmap_crc_cnv.csv")
mut = pd.read_csv(INTERIM_DIR / "depmap_crc_mutation_matrix.csv")

for df in [cohort, gene_effect, expr, cnv, mut]:
    df["ModelID"] = df["ModelID"].astype(str)

# Prefix feature columns
expr_feat = expr.copy()
expr_feat = expr_feat.rename(columns={c: f"expr_{c}" for c in expr_feat.columns if c != "ModelID"})

cnv_feat = cnv.copy()
cnv_feat = cnv_feat.rename(columns={c: f"cnv_{c}" for c in cnv_feat.columns if c != "ModelID"})

mut_feat = mut.copy()

# Feature table without y
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
groups = feature_data["PatientID"].fillna(feature_data["ModelID"]).astype(str)

print(f"Samples: {X_all.shape[0]}")
print(f"Raw features: {X_all.shape[1]}")

cv = GroupKFold(n_splits=5)
k = min(200, X_all.shape[1])

model = Pipeline([
    ("select", SelectKBest(score_func=f_regression, k=k)),
    ("scale", StandardScaler()),
    ("elastic", ElasticNet(alpha=0.05, l1_ratio=0.5, max_iter=10000, random_state=42)),
])

results = []

for gene in TARGET_GENES:
    target_col = find_gene_column(gene_effect.columns, gene)

    if target_col is None:
        print(f"[SKIP] {gene}: target column not found")
        continue

    y_df = gene_effect[["ModelID", target_col]].rename(columns={target_col: "y"})
    data = feature_data[["ModelID", "PatientID", "CellLineName"]].merge(
        y_df, on="ModelID", how="inner"
    )

    y = data["y"].astype(float)

    valid = ~y.isna()
    X = X_all.loc[valid].copy()
    y = y.loc[valid].copy()
    g = groups.loc[valid].copy()

    if len(y) < 20:
        print(f"[SKIP] {gene}: too few samples")
        continue

    pred = cross_val_predict(model, X, y, cv=cv, groups=g)

    mae = mean_absolute_error(y, pred)
    rmse = mean_squared_error(y, pred) ** 0.5
    r2 = r2_score(y, pred)

    baseline_pred = np.repeat(y.mean(), len(y))
    baseline_rmse = mean_squared_error(y, baseline_pred) ** 0.5

    results.append({
        "target_gene": gene,
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

    print(f"\n{gene}")
    print(f"  y range: {y.min():.3f} to {y.max():.3f}")
    print(f"  y std:   {y.std():.3f}")
    print(f"  MAE:     {mae:.4f}")
    print(f"  RMSE:    {rmse:.4f}")
    print(f"  R2:      {r2:.4f}")

results_df = pd.DataFrame(results).sort_values("R2", ascending=False)

out_path = OUT_DIR / "pilot_target_panel_elastic_net_results.csv"
results_df.to_csv(out_path, index=False, encoding="utf-8-sig")

print("\n=== Summary ===")
print(results_df[["target_gene", "y_std", "y_min", "y_max", "RMSE", "baseline_RMSE", "R2"]])

print(f"\nSaved results to: {out_path}")