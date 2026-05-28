from pathlib import Path
import pandas as pd

DATA_DIR = Path("data/raw/depmap")
OUT_DIR = Path("data/interim")
OUT_DIR.mkdir(parents=True, exist_ok=True)

MATRIX_FILES = {
    "gene_effect": "CRISPRGeneEffect.csv",
    "expression": "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
    "cnv": "OmicsCNGene.csv",
}


def read_model_ids(filename):
    path = DATA_DIR / filename
    header = pd.read_csv(path, nrows=0)
    first_col = header.columns[0]

    ids = pd.read_csv(path, usecols=[first_col])
    ids = ids[first_col].astype(str)

    return set(ids)


print("Reading ModelID sets from matrix files...")

id_sets = {}
for key, filename in MATRIX_FILES.items():
    ids = read_model_ids(filename)
    id_sets[key] = ids
    print(f"{key}: {len(ids)} models")

common_ids = set.intersection(*id_sets.values())
print(f"\nCommon models with CRISPR + expression + CNV: {len(common_ids)}")

model = pd.read_csv(DATA_DIR / "Model.csv", low_memory=False)
model["ModelID"] = model["ModelID"].astype(str)

model_common = model[model["ModelID"].isin(common_ids)].copy()
print(f"Models also found in Model.csv: {len(model_common)}")

search_cols = [
    "OncotreeLineage",
    "OncotreePrimaryDisease",
    "OncotreeSubtype",
    "TissueOrigin",
    "SampleCollectionSite",
    "CellLineName",
]

existing_cols = [col for col in search_cols if col in model_common.columns]

combined_text = (
    model_common[existing_cols]
    .fillna("")
    .astype(str)
    .agg(" ".join, axis=1)
    .str.lower()
)

crc_pattern = "colorectal|colon|rectal|rectum|large intestine|bowel"

crc = model_common[combined_text.str.contains(crc_pattern, regex=True)].copy()

print(f"\nCandidate colorectal cancer models: {len(crc)}")

print("\nOncotreeLineage counts:")
print(crc["OncotreeLineage"].value_counts(dropna=False).head(20))

print("\nOncotreePrimaryDisease counts:")
print(crc["OncotreePrimaryDisease"].value_counts(dropna=False).head(30))

keep_cols = [
    "ModelID",
    "PatientID",
    "CellLineName",
    "OncotreeLineage",
    "OncotreePrimaryDisease",
    "OncotreeSubtype",
    "OncotreeCode",
    "PrimaryOrMetastasis",
    "SampleCollectionSite",
    "DepmapModelType",
    "ModelType",
]

keep_cols = [col for col in keep_cols if col in crc.columns]

out_path = OUT_DIR / "depmap_crc_cohort.csv"
crc[keep_cols].to_csv(out_path, index=False, encoding="utf-8-sig")

print(f"\nSaved colorectal cohort to: {out_path}")