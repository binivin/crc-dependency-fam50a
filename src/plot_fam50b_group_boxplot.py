from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

INTERIM_DIR = Path("data/interim")
FIG_DIR = Path("artifacts/figures")
TABLE_DIR = Path("artifacts/tables")
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

TARGET_COL = "FAM50A (9130)"
FEATURE_COL = "FAM50B (26240)"

cohort = pd.read_csv(INTERIM_DIR / "depmap_crc_cohort.csv")
gene_effect = pd.read_csv(INTERIM_DIR / "depmap_crc_gene_effect.csv")
expr = pd.read_csv(INTERIM_DIR / "depmap_crc_expression.csv")

for df in [cohort, gene_effect, expr]:
    df["ModelID"] = df["ModelID"].astype(str)

data = (
    cohort[["ModelID", "CellLineName"]]
    .merge(gene_effect[["ModelID", TARGET_COL]], on="ModelID")
    .merge(expr[["ModelID", FEATURE_COL]], on="ModelID")
)

data = data.rename(columns={
    TARGET_COL: "FAM50A_gene_effect",
    FEATURE_COL: "FAM50B_expression"
}).dropna()

median_expr = data["FAM50B_expression"].median()
data["FAM50B_group"] = data["FAM50B_expression"].apply(
    lambda x: "FAM50B low" if x <= median_expr else "FAM50B high"
)

low = data[data["FAM50B_group"] == "FAM50B low"]["FAM50A_gene_effect"]
high = data[data["FAM50B_group"] == "FAM50B high"]["FAM50A_gene_effect"]

u, p = mannwhitneyu(low, high, alternative="two-sided")

print("FAM50B low group")
print(low.describe())

print("\nFAM50B high group")
print(high.describe())

print(f"\nMann-Whitney U p-value: {p:.6g}")

out_table = TABLE_DIR / "FAM50B_low_high_FAM50A_dependency.csv"
data.to_csv(out_table, index=False, encoding="utf-8-sig")

plt.figure(figsize=(6, 5))
plt.boxplot(
    [low, high],
    labels=["FAM50B low", "FAM50B high"],
    showfliers=True
)

plt.ylabel("FAM50A CRISPR gene effect")
plt.title(f"FAM50A dependency by FAM50B expression group\nMann-Whitney p={p:.2e}")
plt.axhline(-0.5, linestyle="--", linewidth=1)
plt.tight_layout()

out_fig = FIG_DIR / "FAM50A_dependency_by_FAM50B_group.png"
plt.savefig(out_fig, dpi=300)

print(f"\nSaved table to: {out_table}")
print(f"Saved figure to: {out_fig}")