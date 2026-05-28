from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, pearsonr, spearmanr

PROCESSED_DIR = Path("data/processed")
TABLE_DIR = Path("artifacts/tables")
FIG_DIR = Path("artifacts/figures")
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

IN_PATH = PROCESSED_DIR / "tcga_crc_fam50_expression.csv"

print("Reading TCGA FAM50 expression table...")
df = pd.read_csv(IN_PATH)

print(f"Samples: {len(df)}")

# Define FAM50B-low/high by median expression
median_fam50b = df["FAM50B_log2_tpm1"].median()

df["FAM50B_group"] = df["FAM50B_log2_tpm1"].apply(
    lambda x: "FAM50B-low" if x <= median_fam50b else "FAM50B-high"
)

print(f"FAM50B median log2(TPM+1): {median_fam50b:.4f}")

print("\nGroup counts:")
print(df["FAM50B_group"].value_counts())

print("\nProject by FAM50B group:")
project_table = pd.crosstab(df["project_id"], df["FAM50B_group"])
print(project_table)

# FAM50A expression difference by FAM50B group
low = df[df["FAM50B_group"].eq("FAM50B-low")]["FAM50A_log2_tpm1"]
high = df[df["FAM50B_group"].eq("FAM50B-high")]["FAM50A_log2_tpm1"]

u, p = mannwhitneyu(low, high, alternative="two-sided")

print("\nFAM50A expression by FAM50B group:")
print(f"low n = {len(low)}, median = {low.median():.4f}, mean = {low.mean():.4f}")
print(f"high n = {len(high)}, median = {high.median():.4f}, mean = {high.mean():.4f}")
print(f"Mann-Whitney p = {p:.4g}")

# Correlation between FAM50A and FAM50B expression
pearson_r, pearson_p = pearsonr(df["FAM50B_log2_tpm1"], df["FAM50A_log2_tpm1"])
spearman_r, spearman_p = spearmanr(df["FAM50B_log2_tpm1"], df["FAM50A_log2_tpm1"])

print("\nFAM50B vs FAM50A expression correlation:")
print(f"Pearson r = {pearson_r:.4f}, p = {pearson_p:.4g}")
print(f"Spearman r = {spearman_r:.4f}, p = {spearman_p:.4g}")

# Save grouped table
out_table = TABLE_DIR / "tcga_crc_fam50b_groups.csv"
df.to_csv(out_table, index=False, encoding="utf-8-sig")

# Save summary
summary = pd.DataFrame([{
    "n_samples": len(df),
    "fam50b_median_log2_tpm1": median_fam50b,
    "fam50b_low_n": (df["FAM50B_group"] == "FAM50B-low").sum(),
    "fam50b_high_n": (df["FAM50B_group"] == "FAM50B-high").sum(),
    "fam50a_low_group_median": low.median(),
    "fam50a_high_group_median": high.median(),
    "fam50a_group_mannwhitney_p": p,
    "fam50a_fam50b_pearson_r": pearson_r,
    "fam50a_fam50b_pearson_p": pearson_p,
    "fam50a_fam50b_spearman_r": spearman_r,
    "fam50a_fam50b_spearman_p": spearman_p,
}])

summary_path = TABLE_DIR / "tcga_crc_fam50_expression_summary.csv"
summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

# Plot 1: FAM50B distribution
plt.figure(figsize=(7, 5))
plt.hist(df["FAM50B_log2_tpm1"], bins=30)
plt.axvline(median_fam50b, linestyle="--")
plt.xlabel("FAM50B log2(TPM+1)")
plt.ylabel("Number of tumors")
plt.title("TCGA COAD/READ FAM50B expression distribution")
plt.tight_layout()

fig1 = FIG_DIR / "tcga_crc_FAM50B_expression_distribution.png"
plt.savefig(fig1, dpi=300)

# Plot 2: FAM50A by FAM50B group
plt.figure(figsize=(6, 5))
plt.boxplot(
    [low, high],
    tick_labels=["FAM50B-low", "FAM50B-high"],
    showfliers=True
)
plt.ylabel("FAM50A log2(TPM+1)")
plt.title(f"FAM50A expression by FAM50B group\nMann-Whitney p={p:.2e}")
plt.tight_layout()

fig2 = FIG_DIR / "tcga_crc_FAM50A_by_FAM50B_group.png"
plt.savefig(fig2, dpi=300)

# Plot 3: FAM50A vs FAM50B scatter
plt.figure(figsize=(6, 5))
plt.scatter(df["FAM50B_log2_tpm1"], df["FAM50A_log2_tpm1"], alpha=0.7)
plt.xlabel("FAM50B log2(TPM+1)")
plt.ylabel("FAM50A log2(TPM+1)")
plt.title(f"TCGA COAD/READ FAM50A vs FAM50B expression\nPearson r={pearson_r:.3f}")
plt.tight_layout()

fig3 = FIG_DIR / "tcga_crc_FAM50A_vs_FAM50B_expression.png"
plt.savefig(fig3, dpi=300)

print(f"\nSaved grouped table to: {out_table}")
print(f"Saved summary to: {summary_path}")
print(f"Saved figures:")
print(f"  {fig1}")
print(f"  {fig2}")
print(f"  {fig3}")