from pathlib import Path
import pandas as pd
import numpy as np

GDC_DIR = Path("data/raw/gdc/expression")
MANIFEST_PATH = Path("artifacts/manifests/gdc_crc_expression_files.tsv")
OUT_DIR = Path("data/processed")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PATH = OUT_DIR / "tcga_crc_fam50_expression.csv"

TARGET_GENES = ["FAM50A", "FAM50B"]

print("Reading GDC expression metadata...")
meta = pd.read_csv(MANIFEST_PATH, sep="\t")

# Use only primary tumor samples
meta = meta[meta["sample_type"].eq("Primary Tumor")].copy()

print(f"Primary tumor metadata rows: {len(meta)}")

# Map file_name to metadata row
meta_by_filename = meta.set_index("file_name").to_dict(orient="index")

tsv_files = list(GDC_DIR.rglob("*.tsv"))

print(f"Downloaded TSV files found: {len(tsv_files)}")

rows = []

for i, path in enumerate(tsv_files, start=1):
    file_name = path.name

    if file_name not in meta_by_filename:
        # Skip non-primary tumor or unmatched files
        continue

    info = meta_by_filename[file_name]

    df = pd.read_csv(path, sep="\t", comment="#")

    sub = df[df["gene_name"].isin(TARGET_GENES)].copy()

    if sub.empty:
        print(f"[WARNING] No target genes found in {file_name}")
        continue

    values = {
        "file_name": file_name,
        "file_path": str(path),
        "project_id": info.get("project_id"),
        "case_submitter_id": info.get("case_submitter_id"),
        "sample_submitter_id": info.get("sample_submitter_id"),
        "sample_type": info.get("sample_type"),
        "tissue_type": info.get("tissue_type"),
    }

    for gene in TARGET_GENES:
        gene_row = sub[sub["gene_name"].eq(gene)]

        if gene_row.empty:
            values[f"{gene}_tpm"] = np.nan
            values[f"{gene}_log2_tpm1"] = np.nan
            values[f"{gene}_fpkm"] = np.nan
        else:
            r = gene_row.iloc[0]
            tpm = float(r["tpm_unstranded"])
            fpkm = float(r["fpkm_unstranded"])

            values[f"{gene}_tpm"] = tpm
            values[f"{gene}_log2_tpm1"] = np.log2(tpm + 1)
            values[f"{gene}_fpkm"] = fpkm

    rows.append(values)

    if i % 100 == 0:
        print(f"Processed {i} files...")

result = pd.DataFrame(rows)

result = result.sort_values(["project_id", "case_submitter_id", "sample_submitter_id"])

result.to_csv(OUT_PATH, index=False, encoding="utf-8-sig")

print(f"\nSaved FAM50 expression table to: {OUT_PATH}")
print(f"Rows: {len(result)}")

print("\nProject counts:")
print(result["project_id"].value_counts(dropna=False))

print("\nFAM50B log2(TPM+1) summary:")
print(result["FAM50B_log2_tpm1"].describe())

print("\nFAM50A log2(TPM+1) summary:")
print(result["FAM50A_log2_tpm1"].describe())

print("\nFirst rows:")
print(result.head().to_string(index=False))