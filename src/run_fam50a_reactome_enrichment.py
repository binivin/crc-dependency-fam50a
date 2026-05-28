from pathlib import Path
import requests
import pandas as pd

TABLE_DIR = Path("artifacts/tables")
TABLE_DIR.mkdir(parents=True, exist_ok=True)

REACTOME_URL = "https://reactome.org/AnalysisService/identifiers/projection"

GENE_SET_FILES = {
    "all": "FAM50A_pathway_genes_all.txt",
    "stable": "FAM50A_pathway_genes_stable.txt",
    "expression": "FAM50A_pathway_genes_expression.txt",
    "high_confidence": "FAM50A_pathway_genes_high_confidence.txt",
}


def read_genes(path):
    genes = []
    with open(path, "r", encoding="utf-8-sig") as f:
        for line in f:
            gene = line.strip()
            if gene:
                genes.append(gene)
    return genes


def run_reactome(genes, gene_set_name, page_size=100):
    # Reactome recommends a first header line beginning with #
    payload = "#Genes\n" + "\n".join(genes)

    params = {
        "pageSize": page_size,
        "page": 1,
    }

    headers = {
        "Content-Type": "text/plain",
        "Accept": "application/json",
    }

    print(f"\nRunning Reactome enrichment: {gene_set_name}")
    print(f"  n_genes submitted: {len(genes)}")

    response = requests.post(
        REACTOME_URL,
        params=params,
        headers=headers,
        data=payload.encode("utf-8"),
        timeout=120,
    )

    if response.status_code != 200:
        print(response.text[:1000])
        raise RuntimeError(f"Reactome API error: {response.status_code}")

    result = response.json()

    summary = result.get("summary", {})
    token = summary.get("token")
    analysis_type = summary.get("type")

    print(f"  analysis type: {analysis_type}")
    print(f"  token: {token}")

    pathways = result.get("pathways", [])

    rows = []
    for p in pathways:
        entities = p.get("entities", {})
        reactions = p.get("reactions", {})
        species = p.get("species", {})

        rows.append({
            "gene_set": gene_set_name,
            "pathway_name": p.get("name"),
            "stId": p.get("stId"),
            "species": species.get("name") if isinstance(species, dict) else species,
            "entities_found": entities.get("found"),
            "entities_total": entities.get("total"),
            "entities_ratio": entities.get("ratio"),
            "p_value": entities.get("pValue"),
            "fdr": entities.get("fdr"),
            "reactions_found": reactions.get("found"),
            "reactions_total": reactions.get("total"),
            "reactions_ratio": reactions.get("ratio"),
            "token": token,
            "reactome_browser_url": (
                f"https://reactome.org/PathwayBrowser/#/{p.get('stId')}&DTAB=AN&ANALYSIS={token}"
                if token and p.get("stId") else None
            ),
        })

    df = pd.DataFrame(rows)

    if not df.empty:
        df = df.sort_values(["fdr", "p_value"], ascending=True)

    return df, {
        "gene_set": gene_set_name,
        "n_genes_submitted": len(genes),
        "analysis_type": analysis_type,
        "token": token,
        "n_pathways_returned": len(df),
    }


all_results = []
summary_rows = []

for gene_set_name, filename in GENE_SET_FILES.items():
    path = TABLE_DIR / filename

    if not path.exists():
        print(f"[SKIP] Missing file: {path}")
        continue

    genes = read_genes(path)

    if len(genes) < 3:
        print(f"[SKIP] Too few genes for {gene_set_name}: {len(genes)}")
        continue

    df, summary = run_reactome(genes, gene_set_name)

    out_path = TABLE_DIR / f"Reactome_{gene_set_name}_enrichment.csv"
    df.to_csv(out_path, index=False, encoding="utf-8-sig")

    print(f"  saved: {out_path}")

    if not df.empty:
        print("  top pathways:")
        print(
            df[
                ["pathway_name", "entities_found", "entities_total", "p_value", "fdr"]
            ].head(10).to_string(index=False)
        )
    else:
        print("  no pathways returned")

    all_results.append(df)
    summary_rows.append(summary)

if all_results:
    combined = pd.concat(all_results, ignore_index=True)
    combined_out = TABLE_DIR / "Reactome_FAM50A_combined_enrichment.csv"
    combined.to_csv(combined_out, index=False, encoding="utf-8-sig")
    print(f"\nSaved combined Reactome results to: {combined_out}")

summary_df = pd.DataFrame(summary_rows)
summary_out = TABLE_DIR / "Reactome_FAM50A_enrichment_summary.csv"
summary_df.to_csv(summary_out, index=False, encoding="utf-8-sig")

print(f"Saved Reactome summary to: {summary_out}")

print("\nDone.")