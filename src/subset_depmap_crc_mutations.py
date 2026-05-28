from pathlib import Path
import pandas as pd

DATA_DIR = Path("data/raw/depmap")
INTERIM_DIR = Path("data/interim")
INTERIM_DIR.mkdir(parents=True, exist_ok=True)

cohort = pd.read_csv(INTERIM_DIR / "depmap_crc_cohort.csv")
crc_ids = set(cohort["ModelID"].astype(str))

mut_path = DATA_DIR / "OmicsSomaticMutations.csv"
out_path = INTERIM_DIR / "depmap_crc_mutations.csv"

print(f"CRC cohort size: {len(crc_ids)} models")
print("Reading mutation file header...")

header = pd.read_csv(mut_path, nrows=0)
cols = list(header.columns)

print("\nMutation columns:")
print(cols[:50])

possible_id_cols = [
    "ModelID",
    "DepMap_ID",
    "DepMapID",
    "DepmapModelID",
    "DepmapModel",
    "Model",
]

id_col = None
for col in possible_id_cols:
    if col in cols:
        id_col = col
        break

if id_col is None:
    raise ValueError("Could not find a ModelID-like column in OmicsSomaticMutations.csv")

print(f"\nUsing ID column: {id_col}")

chunks = []
total_rows = 0
kept_rows = 0

for chunk in pd.read_csv(mut_path, chunksize=100000, low_memory=False):
    chunk[id_col] = chunk[id_col].astype(str)
    sub = chunk[chunk[id_col].isin(crc_ids)].copy()

    total_rows += len(chunk)
    kept_rows += len(sub)

    if not sub.empty:
        chunks.append(sub)

if chunks:
    result = pd.concat(chunks, axis=0)
    result.to_csv(out_path, index=False, encoding="utf-8-sig")

    print(f"\nTotal mutation rows scanned: {total_rows}")
    print(f"Rows kept: {kept_rows}")
    print(f"Output shape: {result.shape}")
    print(f"Saved to: {out_path}")
else:
    print("\nNo mutation rows matched the CRC cohort.")