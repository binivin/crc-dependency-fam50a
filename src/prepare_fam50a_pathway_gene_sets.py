from pathlib import Path
import re
import pandas as pd

TABLE_DIR = Path("artifacts/tables")
TABLE_DIR.mkdir(parents=True, exist_ok=True)

STABILITY_PATH = TABLE_DIR / "diagnostic_FAM50A_feature_stability.csv"


def extract_source(feature):
    feature = str(feature)

    if feature.startswith("expr_"):
        return "expression"
    if feature.startswith("cnv_"):
        return "cnv"
    if feature.startswith("mut_"):
        return "mutation"

    return "unknown"


def extract_gene_symbol(feature):
    feature = str(feature)

    # expr_FAM50B (26240) -> FAM50B
    # cnv_FAM50B (26240) -> FAM50B
    # mut_FAM50B -> FAM50B
    feature = re.sub(r"^(expr_|cnv_|mut_)", "", feature)
    feature = re.sub(r"\s*\(\d+\)$", "", feature)

    return feature


def is_likely_reactome_mappable_gene(gene):
    gene = str(gene)

    # Very simple filtering to remove obvious non-protein-coding-like symbols.
    # This is not a perfect annotation filter, but it helps avoid clearly unmappable IDs.
    if gene.startswith("RNU"):
        return False
    if gene.startswith("MIR"):
        return False
    if gene.startswith("SNORD"):
        return False
    if gene.startswith("SNORA"):
        return False

    return True


print("Reading FAM50A feature stability table...")
stability = pd.read_csv(STABILITY_PATH)

required_cols = {"feature", "n_folds", "mean_abs_coef", "mean_coef"}
missing = required_cols - set(stability.columns)

if missing:
    raise ValueError(f"Missing columns in stability table: {missing}")

summary = stability.copy()
summary["source"] = summary["feature"].apply(extract_source)
summary["gene_symbol"] = summary["feature"].apply(extract_gene_symbol)
summary["reactome_mappable_candidate"] = summary["gene_symbol"].apply(is_likely_reactome_mappable_gene)

# Put FAM50A and FAM50B at the center of the pathway gene sets
core_genes = pd.DataFrame([
    {
        "feature": "target_FAM50A",
        "n_folds": 5,
        "mean_abs_coef": None,
        "mean_coef": None,
        "source": "target",
        "gene_symbol": "FAM50A",
        "reactome_mappable_candidate": True,
    },
    {
        "feature": "expr_FAM50B (26240)",
        "n_folds": 5,
        "mean_abs_coef": None,
        "mean_coef": None,
        "source": "expression",
        "gene_symbol": "FAM50B",
        "reactome_mappable_candidate": True,
    },
])

summary = pd.concat([core_genes, summary], ignore_index=True)

# Remove duplicate gene/source rows while keeping stronger fold evidence first
summary = summary.sort_values(
    ["n_folds", "mean_abs_coef"],
    ascending=False,
    na_position="last"
)

summary = summary.drop_duplicates(subset=["gene_symbol", "source"], keep="first")

summary_out = TABLE_DIR / "FAM50A_pathway_feature_summary.csv"
summary.to_csv(summary_out, index=False, encoding="utf-8-sig")

print(f"Saved feature summary to: {summary_out}")

# Gene set definitions
all_genes = (
    summary[summary["reactome_mappable_candidate"]]["gene_symbol"]
    .dropna()
    .astype(str)
    .drop_duplicates()
    .sort_values()
    .tolist()
)

stable_genes = (
    summary[
        (summary["n_folds"] >= 2) &
        (summary["reactome_mappable_candidate"])
    ]["gene_symbol"]
    .dropna()
    .astype(str)
    .drop_duplicates()
    .sort_values()
    .tolist()
)

expression_genes = (
    summary[
        (summary["source"].eq("expression")) &
        (summary["reactome_mappable_candidate"])
    ]["gene_symbol"]
    .dropna()
    .astype(str)
    .drop_duplicates()
    .sort_values()
    .tolist()
)

high_confidence_genes = (
    summary[
        (summary["n_folds"] >= 3) &
        (summary["reactome_mappable_candidate"])
    ]["gene_symbol"]
    .dropna()
    .astype(str)
    .drop_duplicates()
    .sort_values()
    .tolist()
)

gene_sets = {
    "FAM50A_pathway_genes_all.txt": all_genes,
    "FAM50A_pathway_genes_stable.txt": stable_genes,
    "FAM50A_pathway_genes_expression.txt": expression_genes,
    "FAM50A_pathway_genes_high_confidence.txt": high_confidence_genes,
}

for filename, genes in gene_sets.items():
    path = TABLE_DIR / filename
    pd.Series(genes).to_csv(path, index=False, header=False)

    print(f"\n{filename}")
    print(f"  n_genes = {len(genes)}")
    print(f"  saved to {path}")
    print(f"  genes: {genes[:30]}")

print("\nDone.")