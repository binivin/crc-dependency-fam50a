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

TARGET_GENE = "KRAS"


def find_gene_column(columns, gene):
    hits = []
    for col in columns:
        col = str(col)
        if col == gene or col.startswith(gene + " ") or col.startswith(gene + " ("):
            hits.append(col)
    if not hits:
        raise ValueError(f"Could not find column for gene: {gene}")
    return hits[0]


print("Reading CRC matrices...")

cohort = pd.read_csv(INTERIM_DIR / "depmap_crc_cohort.csv")
gene_effect = pd.read_csv(INTERIM_DIR / "depmap_crc_gene_effect.csv")
expr = pd.read_csv(INTERIM_DIR / "depmap_crc_expression.csv")
cnv = pd.read_csv(INTERIM_DIR / "depmap_crc_cnv.csv")
mut = pd.read_csv(INTERIM_DIR / "depmap_crc_mutation_matrix.csv")

for df in [cohort, gene_effect, expr, cnv, mut]:
    df["ModelID"] = df["ModelID"].astype(str)

target_col = find_gene_column(gene_effect.columns, TARGET_GENE)

print(f"Target gene: {TARGET_GENE}")
print(f"Target column: {target_col}")

# y: CRISPR dependency score
y_df = gene_effect[["ModelID", target_col]].rename(columns={target_col: "y"})

# Prefix feature columns
expr_feat = expr.copy()
expr_feat = expr_feat.rename(columns={c: f"expr_{c}" for c in expr_feat.columns if c != "ModelID"})

cnv_feat = cnv.copy()
cnv_feat = cnv_feat.rename(columns={c: f"cnv_{c}" for c in cnv_feat.columns if c != "ModelID"})

mut_feat = mut.copy()

# Merge all features
data = (
    y_df
    .merge(cohort[["ModelID", "PatientID", "CellLineName"]], on="ModelID", how="inner")
    .merge(expr_feat, on="ModelID", how="inner")
    .merge(cnv_feat, on="ModelID", how="inner")
    .merge(mut_feat, on="ModelID", how="inner")
)

print(f"Combined data shape: {data.shape}")

# Remove missing target
data = data.dropna(subset=["y"]).copy()
print(f"After removing missing y: {data.shape}")

groups = data["PatientID"].fillna(data["ModelID"]).astype(str)

feature_cols = [
    c for c in data.columns
    if c not in ["ModelID", "PatientID", "CellLineName", "y"]
]

X = data[feature_cols].replace([np.inf, -np.inf], np.nan).fillna(0)
y = data["y"].astype(float)

print(f"Number of samples: {X.shape[0]}")
print(f"Number of raw features: {X.shape[1]}")

# Because n=59 and p is huge, select only top k features inside the pipeline.
k = min(200, X.shape[1])

model = Pipeline([
    ("select", SelectKBest(score_func=f_regression, k=k)),
    ("scale", StandardScaler()),
    ("elastic", ElasticNet(alpha=0.05, l1_ratio=0.5, max_iter=10000, random_state=42)),
])

cv = GroupKFold(n_splits=5)

print("\nRunning 5-fold GroupKFold cross-validation...")
pred = cross_val_predict(model, X, y, cv=cv, groups=groups)

mae = mean_absolute_error(y, pred)
rmse = mean_squared_error(y, pred) ** 0.5
r2 = r2_score(y, pred)

print("\nPilot model performance")
print(f"Target: {TARGET_GENE}")
print(f"MAE:  {mae:.4f}")
print(f"RMSE: {rmse:.4f}")
print(f"R2:   {r2:.4f}")

result = data[["ModelID", "PatientID", "CellLineName", "y"]].copy()
result["pred"] = pred
result["residual"] = result["y"] - result["pred"]

out_path = OUT_DIR / f"pilot_{TARGET_GENE}_elastic_net_predictions.csv"
result.to_csv(out_path, index=False, encoding="utf-8-sig")

print(f"\nSaved predictions to: {out_path}")