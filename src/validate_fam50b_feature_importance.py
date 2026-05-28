from pathlib import Path
import warnings
import numpy as np
import pandas as pd

from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from sklearn.linear_model import ElasticNet, LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore", category=ConvergenceWarning)

INTERIM_DIR = Path("data/interim")
TABLE_DIR = Path("artifacts/tables")
TABLE_DIR.mkdir(parents=True, exist_ok=True)

TARGET_COL = "FAM50A (9130)"
TOP_FEATURE_COL = "FAM50B (26240)"

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

data = (
    cohort[["ModelID", "PatientID", "CellLineName"]]
    .merge(y_df, on="ModelID", how="inner")
    .merge(expr_feat, on="ModelID", how="inner")
    .merge(cnv_feat, on="ModelID", how="inner")
    .merge(mut_feat, on="ModelID", how="inner")
)

data = data.dropna(subset=["y"]).copy()

groups = data["PatientID"].fillna(data["ModelID"]).astype(str)
y = data["y"].astype(float)

fam50b_feature = f"expr_{TOP_FEATURE_COL}"

if fam50b_feature not in data.columns:
    raise ValueError(f"{fam50b_feature} was not found.")

all_feature_cols = [
    c for c in data.columns
    if c not in ["ModelID", "PatientID", "CellLineName", "y"]
]

feature_sets = {
    "FAM50B_expression_only": [fam50b_feature],
    "full_all_features": all_feature_cols,
    "full_without_FAM50B_expression": [c for c in all_feature_cols if c != fam50b_feature],
}

cv = GroupKFold(n_splits=5)
results = []

for label, cols in feature_sets.items():
    print(f"\nEvaluating {label}")
    print(f"  n_features: {len(cols)}")

    X = data[cols].replace([np.inf, -np.inf], np.nan).fillna(0)

    if label == "FAM50B_expression_only":
        model = Pipeline([
            ("scale", StandardScaler()),
            ("linear", LinearRegression()),
        ])
    else:
        k = min(200, X.shape[1])
        model = Pipeline([
            ("var", VarianceThreshold(threshold=0.0)),
            ("select", SelectKBest(score_func=f_regression, k=k)),
            ("scale", StandardScaler()),
            ("elastic", ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=50000, random_state=42)),
        ])

    pred = cross_val_predict(model, X, y, cv=cv, groups=groups)

    mae = mean_absolute_error(y, pred)
    rmse = mean_squared_error(y, pred) ** 0.5
    r2 = r2_score(y, pred)
    corr = np.corrcoef(y, pred)[0, 1]

    results.append({
        "feature_set": label,
        "n_features": X.shape[1],
        "MAE": mae,
        "RMSE": rmse,
        "R2": r2,
        "correlation": corr,
    })

    print(f"  MAE:  {mae:.4f}")
    print(f"  RMSE: {rmse:.4f}")
    print(f"  R2:   {r2:.4f}")
    print(f"  corr: {corr:.4f}")

results_df = pd.DataFrame(results).sort_values("R2", ascending=False)

out_path = TABLE_DIR / "FAM50B_feature_importance_validation.csv"
results_df.to_csv(out_path, index=False, encoding="utf-8-sig")

print("\n=== Feature importance validation ===")
print(results_df.to_string(index=False))

# Median split analysis
fam50b_expr = data[fam50b_feature].astype(float)
median_value = fam50b_expr.median()

low_group = data[fam50b_expr <= median_value].copy()
high_group = data[fam50b_expr > median_value].copy()

low_y = low_group["y"].astype(float)
high_y = high_group["y"].astype(float)

u_stat, p_value = mannwhitneyu(low_y, high_y, alternative="two-sided")

summary = pd.DataFrame([{
    "feature": fam50b_feature,
    "median_expression": median_value,
    "low_group_n": len(low_y),
    "high_group_n": len(high_y),
    "low_group_mean_FAM50A_effect": low_y.mean(),
    "high_group_mean_FAM50A_effect": high_y.mean(),
    "low_group_median_FAM50A_effect": low_y.median(),
    "high_group_median_FAM50A_effect": high_y.median(),
    "mannwhitney_u": u_stat,
    "mannwhitney_p": p_value,
}])

summary_path = TABLE_DIR / "FAM50B_expression_group_comparison.csv"
summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

print("\n=== FAM50B expression group comparison ===")
print(summary.to_string(index=False))

print(f"\nSaved model comparison to: {out_path}")
print(f"Saved group comparison to: {summary_path}")