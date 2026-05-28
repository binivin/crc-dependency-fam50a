from pathlib import Path
import requests
import pandas as pd
import numpy as np

OUT_DIR = Path("data/processed")
TABLE_DIR = Path("artifacts/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

GDC_CASES_URL = "https://api.gdc.cancer.gov/cases"

payload = {
    "filters": {
        "op": "in",
        "content": {
            "field": "project.project_id",
            "value": ["TCGA-COAD", "TCGA-READ"]
        }
    },
    "format": "JSON",
    "fields": ",".join([
        "case_id",
        "submitter_id",
        "project.project_id",
        "demographic.gender",
        "demographic.race",
        "demographic.ethnicity",
        "demographic.vital_status",
        "demographic.days_to_death",
        "diagnoses.age_at_diagnosis",
        "diagnoses.days_to_last_follow_up",
        "diagnoses.primary_diagnosis",
        "diagnoses.tumor_grade",
        "diagnoses.ajcc_pathologic_stage",
        "diagnoses.ajcc_pathologic_t",
        "diagnoses.ajcc_pathologic_n",
        "diagnoses.ajcc_pathologic_m",
        "diagnoses.prior_malignancy",
        "diagnoses.prior_treatment"
    ]),
    "size": "2000"
}

print("Querying GDC cases endpoint for TCGA-COAD/READ clinical metadata...")

response = requests.post(GDC_CASES_URL, json=payload, timeout=120)
response.raise_for_status()

hits = response.json()["data"]["hits"]

print(f"Cases returned: {len(hits)}")

rows = []

for case in hits:
    project = case.get("project", {}) or {}
    demographic = case.get("demographic", {}) or {}
    diagnoses = case.get("diagnoses", []) or [{}]

    # Use first diagnosis record for this first-pass analysis
    diagnosis = diagnoses[0] if diagnoses else {}

    rows.append({
        "case_id": case.get("case_id"),
        "case_submitter_id": case.get("submitter_id"),
        "project_id": project.get("project_id"),
        "gender": demographic.get("gender"),
        "race": demographic.get("race"),
        "ethnicity": demographic.get("ethnicity"),
        "vital_status": demographic.get("vital_status"),
        "days_to_death": demographic.get("days_to_death"),
        "age_at_diagnosis": diagnosis.get("age_at_diagnosis"),
        "days_to_last_follow_up": diagnosis.get("days_to_last_follow_up"),
        "primary_diagnosis": diagnosis.get("primary_diagnosis"),
        "tumor_grade": diagnosis.get("tumor_grade"),
        "ajcc_pathologic_stage": diagnosis.get("ajcc_pathologic_stage"),
        "ajcc_pathologic_t": diagnosis.get("ajcc_pathologic_t"),
        "ajcc_pathologic_n": diagnosis.get("ajcc_pathologic_n"),
        "ajcc_pathologic_m": diagnosis.get("ajcc_pathologic_m"),
        "prior_malignancy": diagnosis.get("prior_malignancy"),
        "prior_treatment": diagnosis.get("prior_treatment"),
    })

clinical = pd.DataFrame(rows)

# Convert survival-like fields
for col in ["days_to_death", "days_to_last_follow_up", "age_at_diagnosis"]:
    clinical[col] = pd.to_numeric(clinical[col], errors="coerce")

clinical["overall_survival_days"] = clinical["days_to_death"].fillna(
    clinical["days_to_last_follow_up"]
)

clinical["death_event"] = clinical["vital_status"].astype(str).str.lower().eq("dead").astype(int)

clinical["age_at_diagnosis_years"] = clinical["age_at_diagnosis"] / 365.25

out_path = OUT_DIR / "tcga_crc_clinical_metadata.csv"
clinical.to_csv(out_path, index=False, encoding="utf-8-sig")

print(f"\nSaved clinical metadata to: {out_path}")

print("\nProject counts:")
print(clinical["project_id"].value_counts(dropna=False))

print("\nVital status counts:")
print(clinical["vital_status"].value_counts(dropna=False))

print("\nAJCC stage counts:")
print(clinical["ajcc_pathologic_stage"].value_counts(dropna=False).head(20))

print("\nOverall survival days summary:")
print(clinical["overall_survival_days"].describe())

print("\nFirst rows:")
print(clinical.head().to_string(index=False))