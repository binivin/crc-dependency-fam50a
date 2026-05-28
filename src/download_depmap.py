import argparse
from pathlib import Path

import pandas as pd
import requests


API = "https://depmap.org/portal/api/download/files"

WANTED = {
    "gene_effect": ["CRISPRGeneEffect.csv"],
    "expr": ["OmicsExpressionProteinCodingGenesTPMLogp1.csv"],
    "cnv": ["OmicsCNGene.csv"],
    "mut": ["OmicsSomaticMutationsMatrixDamaging.csv", "OmicsSomaticMutations.csv"],
    "model": ["Model_v2.csv", "Model.csv"],
    "model_condition": ["ModelCondition.csv"],
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", default="DepMap Public 26Q1")
    parser.add_argument("--out", default="artifacts/manifests/depmap_files.csv")
    args = parser.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("Fetching DepMap file list...")

    headers = {
        "User-Agent": "Mozilla/5.0",
        "Accept": "application/json, text/plain, */*",
        "Referer": "https://depmap.org/portal/data_page/?tab=allData",
    }

    response = requests.get(API, headers=headers, timeout=60)

    if response.status_code == 403:
        print("DepMap API returned 403 Forbidden.")
        print("This means the portal blocked this direct Python request.")
        print("Use manual browser download if this continues.")
        return

    response.raise_for_status()

    files = pd.DataFrame(response.json())

    print("Available columns:")
    print(files.columns.tolist())

    if "releaseName" not in files.columns:
        raise ValueError("releaseName column was not found. Check DepMap API response format.")

    release_files = files[files["releaseName"].eq(args.release)].copy()

    if release_files.empty:
        print(f"No files found for release: {args.release}")
        print("\nAvailable releases:")
        print(files["releaseName"].dropna().unique())
        return

    release_files.to_csv(out_path, index=False, encoding="utf-8-sig")
    print(f"\nSaved full file list to: {out_path}")

    print("\nResolved target files:")
    for key, candidates in WANTED.items():
        hit = release_files[release_files["fileName"].isin(candidates)]

        if hit.empty:
            print(f"{key}: NOT FOUND")
        else:
            row = hit.iloc[0]
            print(f"{key}: {row['fileName']}")

            if "downloadUrl" in row:
                print(f"  URL: {row['downloadUrl']}")


if __name__ == "__main__":
    main()