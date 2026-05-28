from pathlib import Path
import warnings
import numpy as np
import pandas as pd

from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


warnings.filterwarnings("ignore", category=ConvergenceWarning)

INTERIM_DIR = Path("data/interim")
TABLE_DIR = Path("artifacts/tables")
TABLE_DIR.mkdir(parents=True, exist_ok=True)

TARGET_COL = "FAM50A (9130)"
TARGET_GENE = "FAM50A"


def evaluate_model(X, y, groups, label, k=200):
    k = min(k, X.shape[1])

    model = Pipeline([
        ("var", VarianceThreshold(threshold=0.0)),
        ("select", SelectKBest(score_func=f_regression, k=k)),
        ("scale", StandardScaler()),
        ("elastic", ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=50000, random_state=42)),
    ])

    cv = GroupKFold(n_splits=5)

    pred = cross_val_predict(model, X, y, cv=cv, groups=groups)

    mae = mean_absolute_error(y, pred)
    rmse = mean_squared_error(y, pred) ** 0.5
    r2 = r2_score(y, pred)
    corr = np.corrcoef(y, pred)[0, 1]

    baseline_pred = np.repeat(y.mean(), len(y))
    baseline_rmse = mean_squared_error(y, baseline_pred) ** 0.5

    return {
        "feature_set": label,
        "n_samples": X.shape[0],
        "n_features": X.shape[1],
        "MAE": mae,
        "RMSE": rmse,
        "baseline_RMSE": baseline_rmse,
        "R2": r2,
        "correlation": corr,
    }


print("Reading data...")

cohort = pd.read_csv(INTERIM_DIR / "depmap_crc_cohort.csv")
gene_effect = pd.read_csv(INTERIM_DIR / "depmap_crc_gene_effect.csv")
expr = pd.read_csv(INTERIM_DIR / "depmap_crc_expression.csv")
cnv = pd.read_csv(INTERIM_DIR / "depmap_crc_cnv.csv")
mut = pd.read_csv(INTERIM_DIR / "depmap_crc_mutation_matrix.csv")

for df in [cohort, gene_effect, expr, cnv, mut]:
    df["ModelID"] = df["ModelID"].astype(str)

y_df = gene_effect[["ModelID", TARGET_COL]].rename(columns={TARGET_COL: "y"})

expr_feat = expr.rename(columns={c: f"expr_{c}" for c in expr.columns if c != "ModelID"})
cnv_feat = cnv.rename(columns={c: f"cnv_{c}" for c in cnv.columns if c != "ModelID"})
mut_feat = mut.copy()

base = (
    cohort[["ModelID", "PatientID", "CellLineName"]]
    .merge(y_df, on="ModelID", how="inner")
)

data_expr = base.merge(expr_feat, on="ModelID", how="inner")
data_cnv = base.merge(cnv_feat, on="ModelID", how="inner")
data_mut = base.merge(mut_feat, on="ModelID", how="inner")
data_expr_cnv = base.merge(expr_feat, on="ModelID", how="inner").merge(cnv_feat, on="ModelID", how="inner")
data_all = data_expr_cnv.merge(mut_feat, on="ModelID", how="inner")

datasets = {
    "expression_only": data_expr,
    "cnv_only": data_cnv,
    "mutation_only": data_mut,
    "expression_cnv": data_expr_cnv,
    "expression_cnv_mutation": data_all,
}

results = []

for label, data in datasets.items():
    data = data.dropna(subset=["y"]).copy()

    feature_cols = [
        c for c in data.columns
        if c not in ["ModelID", "PatientID", "CellLineName", "y"]
    ]

    X = data[feature_cols].replace([np.inf, -np.inf], np.nan).fillna(0)
    y = data["y"].astype(float)
    groups = data["PatientID"].fillna(data["ModelID"]).astype(str)

    print(f"\nEvaluating {label}")
    print(f"  samples: {X.shape[0]}")
    print(f"  features: {X.shape[1]}")

    res = evaluate_model(X, y, groups, label)
    results.append(res)

    print(f"  R2: {res['R2']:.4f}")
    print(f"  corr: {res['correlation']:.4f}")
    print(f"  RMSE: {res['RMSE']:.4f}")

results_df = pd.DataFrame(results).sort_values("R2", ascending=False)

out_path = TABLE_DIR / "FAM50A_ablation_results.csv"
results_df.to_csv(out_path, index=False, encoding="utf-8-sig")

print("\n=== Ablation summary ===")
print(results_df.to_string(index=False))

print(f"\nSaved to: {out_path}")