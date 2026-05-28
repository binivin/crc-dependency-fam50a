from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error, r2_score


warnings.filterwarnings("ignore", category=ConvergenceWarning)

INTERIM_DIR = Path("data/interim")
TABLE_DIR = Path("artifacts/tables")
FIG_DIR = Path("artifacts/figures")
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

TARGET_COL = "FAM50A (9130)"
TARGET_GENE = "FAM50A"
N_PERMUTATIONS = 50
RANDOM_SEED = 42


def build_model(n_features):
    k = min(200, n_features)

    return Pipeline([
        ("var", VarianceThreshold(threshold=0.0)),
        ("select", SelectKBest(score_func=f_regression, k=k)),
        ("scale", StandardScaler()),
        ("elastic", ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=50000, random_state=42)),
    ])


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

feature_cols = [
    c for c in data.columns
    if c not in ["ModelID", "PatientID", "CellLineName", "y"]
]

X = data[feature_cols].replace([np.inf, -np.inf], np.nan).fillna(0)
y = data["y"].astype(float).reset_index(drop=True)
groups = data["PatientID"].fillna(data["ModelID"]).astype(str).reset_index(drop=True)

print(f"Samples: {X.shape[0]}")
print(f"Features: {X.shape[1]}")

cv = GroupKFold(n_splits=5)

# Observed model
print("\nRunning observed model...")
model = build_model(X.shape[1])
observed_pred = cross_val_predict(model, X, y, cv=cv, groups=groups)

observed_rmse = mean_squared_error(y, observed_pred) ** 0.5
observed_r2 = r2_score(y, observed_pred)

print(f"Observed RMSE: {observed_rmse:.4f}")
print(f"Observed R2:   {observed_r2:.4f}")

# Permutation test
rng = np.random.default_rng(RANDOM_SEED)
perm_results = []

print(f"\nRunning {N_PERMUTATIONS} label permutations...")

for i in range(1, N_PERMUTATIONS + 1):
    y_perm = pd.Series(rng.permutation(y.values))

    model = build_model(X.shape[1])
    perm_pred = cross_val_predict(model, X, y_perm, cv=cv, groups=groups)

    perm_rmse = mean_squared_error(y_perm, perm_pred) ** 0.5
    perm_r2 = r2_score(y_perm, perm_pred)

    perm_results.append({
        "iteration": i,
        "RMSE": perm_rmse,
        "R2": perm_r2,
    })

    print(f"[{i:02d}/{N_PERMUTATIONS}] permuted R2 = {perm_r2:.4f}")

perm_df = pd.DataFrame(perm_results)

# Empirical p-value
p_value = (1 + (perm_df["R2"] >= observed_r2).sum()) / (N_PERMUTATIONS + 1)

summary = pd.DataFrame([{
    "target_gene": TARGET_GENE,
    "observed_R2": observed_r2,
    "observed_RMSE": observed_rmse,
    "permutation_mean_R2": perm_df["R2"].mean(),
    "permutation_max_R2": perm_df["R2"].max(),
    "permutation_min_R2": perm_df["R2"].min(),
    "n_permutations": N_PERMUTATIONS,
    "empirical_p_value": p_value,
}])

perm_path = TABLE_DIR / "FAM50A_permutation_results.csv"
summary_path = TABLE_DIR / "FAM50A_permutation_summary.csv"

perm_df.to_csv(perm_path, index=False, encoding="utf-8-sig")
summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

print("\n=== Permutation summary ===")
print(summary.to_string(index=False))

# Plot permutation distribution
plt.figure(figsize=(7, 5))
plt.hist(perm_df["R2"], bins=15)
plt.axvline(observed_r2, linestyle="--", label=f"Observed R2 = {observed_r2:.3f}")
plt.xlabel("R2 under permuted labels")
plt.ylabel("Count")
plt.title("FAM50A permutation test")
plt.legend()
plt.tight_layout()

fig_path = FIG_DIR / "FAM50A_permutation_test.png"
plt.savefig(fig_path, dpi=300)

print(f"\nSaved permutation results to: {perm_path}")
print(f"Saved summary to: {summary_path}")
print(f"Saved figure to: {fig_path}")