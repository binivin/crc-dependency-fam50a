from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

INTERIM_DIR = Path("data/interim")
TABLE_DIR = Path("artifacts/tables")
FIG_DIR = Path("artifacts/figures")
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

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

data = (
    cohort[["ModelID", "PatientID", "CellLineName"]]
    .merge(gene_effect[["ModelID", TARGET_COL]], on="ModelID", how="inner")
    .merge(expr[["ModelID", TOP_FEATURE_COL]], on="ModelID", how="inner")
)

data = data.rename(columns={
    TARGET_COL: "FAM50A_dependency",
    TOP_FEATURE_COL: "FAM50B_expression"
})

data = data.dropna().copy()

x = data["FAM50B_expression"].astype(float)
y = data["FAM50A_dependency"].astype(float)

pearson_r, pearson_p = pearsonr(x, y)
spearman_r, spearman_p = spearmanr(x, y)

print("\n=== FAM50B expression vs FAM50A dependency ===")
print(f"n = {len(data)}")
print(f"Pearson r  = {pearson_r:.4f}, p = {pearson_p:.4g}")
print(f"Spearman r = {spearman_r:.4f}, p = {spearman_p:.4g}")

# Save table
out_table = TABLE_DIR / "FAM50A_FAM50B_relationship.csv"
data.to_csv(out_table, index=False, encoding="utf-8-sig")

# Plot
plt.figure(figsize=(6, 5))
plt.scatter(x, y)

# regression line
coef = np.polyfit(x, y, deg=1)
line_x = np.linspace(x.min(), x.max(), 100)
line_y = coef[0] * line_x + coef[1]
plt.plot(line_x, line_y, linestyle="--")

plt.xlabel("FAM50B expression")
plt.ylabel("FAM50A CRISPR gene effect")
plt.title(f"FAM50B expression vs FAM50A dependency\nPearson r={pearson_r:.3f}, p={pearson_p:.3g}")
plt.tight_layout()

out_fig = FIG_DIR / "FAM50A_FAM50B_relationship.png"
plt.savefig(out_fig, dpi=300)

print(f"\nSaved table to: {out_table}")
print(f"Saved figure to: {out_fig}")

# Also inspect whether FAM50A or FAM50B mutation exists in these cell lines
for gene in ["FAM50A", "FAM50B"]:
    col = f"mut_{gene}"
    if col in mut.columns:
        tmp = cohort[["ModelID", "CellLineName"]].merge(mut[["ModelID", col]], on="ModelID", how="inner")
        print(f"\n{col} mutated models:")
        print(tmp[tmp[col] == 1][["ModelID", "CellLineName"]].to_string(index=False))
    else:
        print(f"\n{col} not found in mutation matrix.")