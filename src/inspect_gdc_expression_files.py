from pathlib import Path
import pandas as pd

GDC_DIR = Path("data/raw/gdc/expression")

print("Searching for downloaded GDC expression TSV files...")

tsv_files = list(GDC_DIR.rglob("*.tsv"))

print(f"Number of TSV files found: {len(tsv_files)}")

if not tsv_files:
    raise FileNotFoundError("No .tsv files found under data/raw/gdc/expression")

for path in tsv_files[:5]:
    print("\n" + "=" * 80)
    print(f"File: {path}")

    # GDC STAR-Counts files usually have comment/header lines.
    df = pd.read_csv(path, sep="\t", comment="#", nrows=10)

    print("\nColumns:")
    print(df.columns.tolist())

    print("\nFirst rows:")
    print(df.head())

    full = pd.read_csv(path, sep="\t", comment="#")

    if "gene_name" in full.columns:
        target = full[full["gene_name"].isin(["FAM50A", "FAM50B"])].copy()

        print("\nFAM50A/FAM50B rows:")
        if target.empty:
            print("No FAM50A/FAM50B rows found.")
        else:
            print(target.to_string(index=False))
    else:
        print("\nNo gene_name column found.")