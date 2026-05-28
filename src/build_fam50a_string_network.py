from pathlib import Path
import re
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

TABLE_DIR = Path("artifacts/tables")
FIG_DIR = Path("artifacts/figures")
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

STABILITY_PATH = TABLE_DIR / "diagnostic_FAM50A_feature_stability.csv"


def extract_gene_symbol(feature):
    feature = str(feature)
    feature = re.sub(r"^(expr_|cnv_|mut_)", "", feature)
    feature = re.sub(r"\s*\(\d+\)$", "", feature)
    return feature


print("Reading stable features...")
stability = pd.read_csv(STABILITY_PATH)

# 반복적으로 선택된 feature만 사용
stable = stability[stability["n_folds"] >= 2].copy()

genes = ["FAM50A", "FAM50B"]

for feat in stable["feature"].head(30):
    gene = extract_gene_symbol(feat)
    if gene not in genes:
        genes.append(gene)

print(f"Gene list size: {len(genes)}")
print(genes)

gene_list_path = TABLE_DIR / "FAM50A_STRING_input_genes.txt"
pd.Series(genes).to_csv(gene_list_path, index=False, header=False)

print(f"\nSaved input genes to: {gene_list_path}")

# STRING API
url = "https://string-db.org/api/tsv/network"

params = {
    "identifiers": "%0d".join(genes),
    "species": 9606,
    "required_score": 400,
    "network_type": "functional",
    "caller_identity": "crc_dependency_project"
}

print("\nQuerying STRING API...")
response = requests.post(url, data=params, timeout=60)

if response.status_code != 200:
    raise RuntimeError(f"STRING API error: {response.status_code}\n{response.text[:500]}")

text = response.text.strip()

if not text:
    print("No STRING network returned.")
    raise SystemExit

rows = [line.split("\t") for line in text.split("\n")]
cols = rows[0]
df = pd.DataFrame(rows[1:], columns=cols)

out_path = TABLE_DIR / "FAM50A_STRING_network.tsv"
df.to_csv(out_path, sep="\t", index=False, encoding="utf-8-sig")

print(f"Saved STRING network to: {out_path}")
print(f"Edges returned: {len(df)}")

score_col = "score" if "score" in df.columns else "combined_score"

G = nx.Graph()

for _, row in df.iterrows():
    a = row["preferredName_A"]
    b = row["preferredName_B"]
    score = float(row[score_col])
    G.add_edge(a, b, weight=score)

print(f"\nNetwork nodes: {G.number_of_nodes()}")
print(f"Network edges: {G.number_of_edges()}")

degree = dict(G.degree())
betweenness = nx.betweenness_centrality(G, weight="weight")

centrality = pd.DataFrame({
    "gene": list(G.nodes()),
    "degree": [degree[g] for g in G.nodes()],
    "betweenness": [betweenness[g] for g in G.nodes()],
})

centrality = centrality.sort_values(["degree", "betweenness"], ascending=False)

centrality_path = TABLE_DIR / "FAM50A_STRING_centrality.csv"
centrality.to_csv(centrality_path, index=False, encoding="utf-8-sig")

print(f"Saved centrality table to: {centrality_path}")

print("\nTop central genes:")
print(centrality.head(20).to_string(index=False))

plt.figure(figsize=(9, 7))

pos = nx.spring_layout(G, seed=42, weight="weight")

node_sizes = []
for node in G.nodes():
    if node in ["FAM50A", "FAM50B"]:
        node_sizes.append(900)
    else:
        node_sizes.append(300 + 80 * degree[node])

nx.draw_networkx_edges(G, pos, alpha=0.4)
nx.draw_networkx_nodes(G, pos, node_size=node_sizes)
nx.draw_networkx_labels(G, pos, font_size=8)

plt.title("STRING functional network around FAM50A model features")
plt.axis("off")
plt.tight_layout()

fig_path = FIG_DIR / "FAM50A_STRING_network.png"
plt.savefig(fig_path, dpi=300)

print(f"Saved network figure to: {fig_path}")