from pathlib import Path
import requests
import pandas as pd

OUT_DIR = Path("artifacts/manifests")
OUT_DIR.mkdir(parents=True, exist_ok=True)

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"

payload = {
    "filters": {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": ["TCGA-COAD", "TCGA-READ"]
                }
            },
            {
                "op": "=",
                "content": {
                    "field": "files.data_type",
                    "value": "Gene Expression Quantification"
                }
            },
            {
                "op": "=",
                "content": {
                    "field": "analysis.workflow_type",
                    "value": "STAR - Counts"
                }
            }
        ]
    },
    "format": "JSON",
    "fields": ",".join([
        "file_id",
        "file_name",
        "file_size",
        "md5sum",
        "data_category",
        "data_type",
        "analysis.workflow_type",
        "cases.submitter_id",
        "cases.case_id",
        "cases.project.project_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_id",
        "cases.samples.sample_type",
        "cases.samples.tissue_type",
        "cases.samples.tumor_descriptor"
    ]),
    "size": "5000"
}

print("Querying GDC files endpoint for TCGA-COAD/READ STAR-Counts expression files...")

response = requests.post(GDC_FILES_URL, json=payload, timeout=120)
response.raise_for_status()

data = response.json()["data"]["hits"]

print(f"Number of files returned: {len(data)}")

rows = []

for item in data:
    cases = item.get("cases", [])

    if not cases:
        rows.append({
            "file_id": item.get("file_id"),
            "file_name": item.get("file_name"),
            "file_size": item.get("file_size"),
            "md5sum": item.get("md5sum"),
            "project_id": None,
            "case_submitter_id": None,
            "sample_submitter_id": None,
            "sample_type": None,
            "tissue_type": None,
        })
        continue

    for case in cases:
        samples = case.get("samples", [])

        if not samples:
            rows.append({
                "file_id": item.get("file_id"),
                "file_name": item.get("file_name"),
                "file_size": item.get("file_size"),
                "md5sum": item.get("md5sum"),
                "project_id": case.get("project", {}).get("project_id"),
                "case_submitter_id": case.get("submitter_id"),
                "sample_submitter_id": None,
                "sample_type": None,
                "tissue_type": None,
            })
            continue

        for sample in samples:
            rows.append({
                "file_id": item.get("file_id"),
                "file_name": item.get("file_name"),
                "file_size": item.get("file_size"),
                "md5sum": item.get("md5sum"),
                "project_id": case.get("project", {}).get("project_id"),
                "case_submitter_id": case.get("submitter_id"),
                "sample_submitter_id": sample.get("submitter_id"),
                "sample_type": sample.get("sample_type"),
                "tissue_type": sample.get("tissue_type"),
                "tumor_descriptor": sample.get("tumor_descriptor"),
                "workflow_type": item.get("analysis", {}).get("workflow_type"),
                "data_type": item.get("data_type"),
            })

df = pd.DataFrame(rows)

out_tsv = OUT_DIR / "gdc_crc_expression_files.tsv"
df.to_csv(out_tsv, sep="\t", index=False, encoding="utf-8-sig")

print(f"Saved metadata table to: {out_tsv}")

print("\nProject counts:")
print(df["project_id"].value_counts(dropna=False))

print("\nSample type counts:")
print(df["sample_type"].value_counts(dropna=False))

# For first patient-level analysis, focus on Primary Tumor samples.
primary = df[df["sample_type"].eq("Primary Tumor")].copy()

print(f"\nPrimary tumor files: {len(primary)}")

manifest = primary[["file_id", "file_name", "md5sum", "file_size"]].copy()
manifest = manifest.rename(columns={
    "file_id": "id",
    "file_name": "filename",
    "file_size": "size"
})

manifest_out = OUT_DIR / "gdc_crc_expression_manifest.txt"
manifest.to_csv(manifest_out, sep="\t", index=False)

print(f"Saved GDC download manifest to: {manifest_out}")

print("\nDone.")