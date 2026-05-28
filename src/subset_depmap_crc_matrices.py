from pathlib import Path
import pandas as pd

DATA_DIR = Path("data/raw/depmap")
INTERIM_DIR = Path("data/interim")
INTERIM_DIR.mkdir(parents=True, exist_ok=True)

cohort = pd.read_csv(INTERIM_DIR / "depmap_crc_cohort.csv")
crc_ids = set(cohort["ModelID"].astype(str))

print(f"CRC cohort size: {len(crc_ids)} models")

FILES = {
    "gene_effect": "CRISPRGeneEffect.csv",
    "expression": "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
    "cnv": "OmicsCNGene.csv",
}

OUT_FILES = {
    "gene_effect": "depmap_crc_gene_effect.csv",
    "expression": "depmap_crc_expression.csv",
    "cnv": "depmap_crc_cnv.csv",
}


def subset_matrix(input_file, output_file, chunksize=100):
    in_path = DATA_DIR / input_file
    out_path = INTERIM_DIR / output_file

    print(f"\nSubsetting {input_file}...")

    chunks = []
    total_rows = 0
    kept_rows = 0

    for chunk in pd.read_csv(in_path, chunksize=chunksize, low_memory=False):
        first_col = chunk.columns[0]
        chunk[first_col] = chunk[first_col].astype(str)

        sub = chunk[chunk[first_col].isin(crc_ids)].copy()

        total_rows += len(chunk)
        kept_rows += len(sub)

        if not sub.empty:
            chunks.append(sub)

    if not chunks:
        print(f"No matching rows found for {input_file}")
        return

    result = pd.concat(chunks, axis=0)

    # Rename first column to ModelID for consistency
    first_col = result.columns[0]
    result = result.rename(columns={first_col: "ModelID"})

    result.to_csv(out_path, index=False, encoding="utf-8-sig")

    print(f"Total rows scanned: {total_rows}")
    print(f"Rows kept: {kept_rows}")
    print(f"Output shape: {result.shape}")
    print(f"Saved to: {out_path}")


for key, input_file in FILES.items():
    subset_matrix(input_file, OUT_FILES[key])

print("\nDone.")