from pathlib import Path
import pandas as pd

DATA_DIR = Path("data/raw/depmap")

files = {
    "gene_effect": "CRISPRGeneEffect.csv",
    "expression": "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
    "cnv": "OmicsCNGene.csv",
    "mutation": "OmicsSomaticMutations.csv",
    "model": "Model.csv",
    "model_condition": "ModelCondition.csv",
}

print("Checking files...\n")

for key, fname in files.items():
    path = DATA_DIR / fname
    if path.exists():
        size_mb = path.stat().st_size / 1024 / 1024
        print(f"[OK] {key}: {fname} ({size_mb:.1f} MB)")
    else:
        print(f"[MISSING] {key}: {fname}")

print("\nInspecting table shapes and ID columns...\n")

# Large matrix files: read only header and first column
matrix_files = ["CRISPRGeneEffect.csv", "OmicsExpressionProteinCodingGenesTPMLogp1.csv", "OmicsCNGene.csv"]

model_ids = {}

for fname in matrix_files:
    path = DATA_DIR / fname

    if not path.exists():
        continue

    header = pd.read_csv(path, nrows=0)
    first_col = header.columns[0]

    ids = pd.read_csv(path, usecols=[0])
    ids_set = set(ids[first_col].astype(str))

    model_ids[fname] = ids_set

    print(f"{fname}")
    print(f"  first column: {first_col}")
    print(f"  number of rows/models: {len(ids_set)}")
    print(f"  number of columns/features: {len(header.columns) - 1}")
    print()

# Model metadata
model_path = DATA_DIR / "Model.csv"
if model_path.exists():
    model = pd.read_csv(model_path, low_memory=False)
    print("Model.csv columns:")
    print(model.columns.tolist()[:30])
    print(f"Model.csv rows: {len(model)}")

    possible_id_cols = ["ModelID", "ModelID"]
    id_col = None
    for col in possible_id_cols:
        if col in model.columns:
            id_col = col
            break

    if id_col is not None:
        model_set = set(model[id_col].astype(str))
        print(f"\nUsing metadata ID column: {id_col}")
        for fname, ids_set in model_ids.items():
            overlap = len(ids_set & model_set)
            print(f"Overlap between {fname} and Model.csv: {overlap}")
    else:
        print("\nCould not find ModelID column in Model.csv.")

print("\nPairwise overlap among matrix files:")
names = list(model_ids.keys())
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = names[i], names[j]
        overlap = len(model_ids[a] & model_ids[b])
        print(f"{a} ∩ {b}: {overlap}")