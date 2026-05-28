from pathlib import Path
import pandas as pd

INTERIM_DIR = Path("data/interim")

mut_path = INTERIM_DIR / "depmap_crc_mutations.csv"
cohort_path = INTERIM_DIR / "depmap_crc_cohort.csv"
out_path = INTERIM_DIR / "depmap_crc_mutation_matrix.csv"

print("Reading mutation table...")
mut = pd.read_csv(mut_path, low_memory=False)

print("Reading CRC cohort...")
cohort = pd.read_csv(cohort_path)
crc_ids = cohort["ModelID"].astype(str).tolist()

mut["ModelID"] = mut["ModelID"].astype(str)

# Keep rows with valid gene symbols
mut = mut.dropna(subset=["HugoSymbol"])
mut["HugoSymbol"] = mut["HugoSymbol"].astype(str)

print(f"Mutation rows after removing missing HugoSymbol: {len(mut)}")

# Optional: focus on likely functional mutations
impact_cols = ["VepImpact", "MolecularConsequence"]

if "VepImpact" in mut.columns:
    print("\nVepImpact counts:")
    print(mut["VepImpact"].value_counts(dropna=False).head(20))

# For MVP, use all observed non-silent gene-level variants.
# Later we can restrict this to HIGH/MODERATE or damaging mutations only.
gene_mut = (
    mut[["ModelID", "HugoSymbol"]]
    .drop_duplicates()
    .assign(value=1)
)

matrix = gene_mut.pivot_table(
    index="ModelID",
    columns="HugoSymbol",
    values="value",
    fill_value=0,
    aggfunc="max"
)

# Ensure all 59 CRC models are present
matrix = matrix.reindex(crc_ids).fillna(0).astype(int)

# Rename mutation columns
matrix.columns = [f"mut_{gene}" for gene in matrix.columns]

matrix = matrix.reset_index()

matrix.to_csv(out_path, index=False, encoding="utf-8-sig")

print(f"\nOutput shape: {matrix.shape}")
print(f"Saved to: {out_path}")

print("\nExample columns:")
print(matrix.columns[:20].tolist())

print("\nMutation burden per model:")
burden = matrix.drop(columns=["ModelID"]).sum(axis=1)
print(burden.describe())