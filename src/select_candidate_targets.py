from pathlib import Path
import pandas as pd
import numpy as np
import re

INTERIM_DIR = Path("data/interim")
OUT_DIR = Path("artifacts/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)

gene_effect_path = INTERIM_DIR / "depmap_crc_gene_effect.csv"
out_path = OUT_DIR / "candidate_dependency_targets.csv"

print("Reading CRC gene effect matrix...")
ge = pd.read_csv(gene_effect_path)

ge["ModelID"] = ge["ModelID"].astype(str)

score_cols = [c for c in ge.columns if c != "ModelID"]

scores = ge[score_cols].apply(pd.to_numeric, errors="coerce")

def extract_gene_symbol(col):
    col = str(col)
    return re.sub(r"\s*\(\d+\)$", "", col)

stats = pd.DataFrame({
    "target_column": score_cols,
    "gene_symbol": [extract_gene_symbol(c) for c in score_cols],
    "y_mean": scores.mean(axis=0).values,
    "y_std": scores.std(axis=0).values,
    "y_min": scores.min(axis=0).values,
    "y_max": scores.max(axis=0).values,
    "n_missing": scores.isna().sum(axis=0).values,
    "n_dep_lt_m05": (scores < -0.5).sum(axis=0).values,
    "n_strong_lt_m1": (scores < -1.0).sum(axis=0).values,
    "n_neutral_gt_m02": (scores > -0.2).sum(axis=0).values,
    "frac_dep_lt_m05": (scores < -0.5).mean(axis=0).values,
    "frac_strong_lt_m1": (scores < -1.0).mean(axis=0).values,
})

# Candidate rules:
# - variation should be large
# - at least some sensitive cell lines
# - at least some neutral/resistant cell lines
# - avoid genes that are strongly essential in almost every cell line
candidates = stats[
    (stats["n_missing"] == 0) &
    (stats["y_std"] >= 0.30) &
    (stats["y_min"] <= -0.8) &
    (stats["n_dep_lt_m05"] >= 5) &
    (stats["n_neutral_gt_m02"] >= 5) &
    (stats["frac_strong_lt_m1"] <= 0.70)
].copy()

# Higher score = more useful for pilot modeling
candidates["selection_score"] = (
    candidates["y_std"]
    + 0.2 * candidates["n_dep_lt_m05"]
    + 0.1 * candidates["n_neutral_gt_m02"]
)

candidates = candidates.sort_values(
    ["selection_score", "y_std", "n_dep_lt_m05"],
    ascending=False
)

candidates.to_csv(out_path, index=False, encoding="utf-8-sig")

print(f"\nTotal genes checked: {len(stats)}")
print(f"Candidate targets found: {len(candidates)}")
print(f"Saved to: {out_path}")

print("\nTop 30 candidate targets:")
print(
    candidates[
        [
            "gene_symbol",
            "target_column",
            "y_std",
            "y_min",
            "y_max",
            "n_dep_lt_m05",
            "n_strong_lt_m1",
            "n_neutral_gt_m02",
            "frac_strong_lt_m1",
        ]
    ].head(30).to_string(index=False)
)