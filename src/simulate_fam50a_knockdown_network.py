from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

PROCESSED_DIR = Path("data/processed")
TABLE_DIR = Path("artifacts/tables")
FIG_DIR = Path("artifacts/figures")

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

VULN_PATH = PROCESSED_DIR / "tcga_crc_fam50a_predicted_vulnerability.csv"
STRING_PATH = TABLE_DIR / "FAM50A_STRING_network.tsv"

print("Reading TCGA predicted vulnerability table...")
tcga = pd.read_csv(VULN_PATH)

print(f"TCGA rows: {len(tcga)}")

required_cols = [
    "case_submitter_id",
    "project_id",
    "FAM50B_log2_tpm1",
    "predicted_FAM50A_gene_effect",
    "predicted_FAM50A_vulnerability_score",
    "FAM50A_predicted_vulnerability_group",
]

missing = [c for c in required_cols if c not in tcga.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Case-level table
case_df = (
    tcga
    .groupby("case_submitter_id", as_index=False)
    .agg(
        project_id=("project_id", "first"),
        FAM50B_log2_tpm1=("FAM50B_log2_tpm1", "median"),
        predicted_FAM50A_gene_effect=("predicted_FAM50A_gene_effect", "median"),
        predicted_FAM50A_vulnerability_score=("predicted_FAM50A_vulnerability_score", "median"),
        FAM50A_predicted_vulnerability_group=("FAM50A_predicted_vulnerability_group", "first"),
    )
)

print(f"Case-level rows: {len(case_df)}")

# Use only positive vulnerability values for fitness-loss simulation.
# Larger score means stronger estimated vulnerability.
case_df["vulnerability_positive"] = case_df["predicted_FAM50A_vulnerability_score"].clip(lower=0)

print("\nVulnerability score summary:")
print(case_df["vulnerability_positive"].describe())

# ---------------------------------------------------------------------
# Build STRING network
# ---------------------------------------------------------------------

print("\nReading STRING network...")

G = nx.Graph()

if STRING_PATH.exists():
    string_df = pd.read_csv(STRING_PATH, sep="\t")

    if len(string_df) > 0:
        score_col = "score" if "score" in string_df.columns else "combined_score"

        for _, row in string_df.iterrows():
            a = row["preferredName_A"]
            b = row["preferredName_B"]
            score = float(row[score_col])

            # Normalize STRING score if needed.
            # STRING API score is usually 0~1 in this output.
            if score > 1:
                score = score / 1000.0

            G.add_edge(a, b, weight=score)

# Safety fallback if STRING graph is too small or missing
if G.number_of_nodes() == 0:
    G.add_edge("FAM50A", "FAM50B", weight=1.0)

print(f"Network nodes: {G.number_of_nodes()}")
print(f"Network edges: {G.number_of_edges()}")
print(f"Nodes: {list(G.nodes())}")
print(f"Edges: {list(G.edges(data=True))}")

# ---------------------------------------------------------------------
# Knockdown propagation model
# ---------------------------------------------------------------------

def propagate_knockdown(G, source="FAM50A", alpha=1.0, beta=0.30, n_steps=5):
    """
    Simple network propagation model.

    loss[node] ranges from 0 to 1.
    alpha: direct knockdown intensity on source node.
    beta: propagation strength through network edges.
    n_steps: number of propagation iterations.

    This is a toy simulation, not a mechanistic biochemical model.
    """
    nodes = list(G.nodes())
    idx = {node: i for i, node in enumerate(nodes)}

    W = np.zeros((len(nodes), len(nodes)))

    for a, b, data in G.edges(data=True):
        weight = float(data.get("weight", 1.0))
        W[idx[a], idx[b]] = weight
        W[idx[b], idx[a]] = weight

    # Row-normalize weights to avoid uncontrolled growth
    row_sums = W.sum(axis=1, keepdims=True)
    W_norm = np.divide(W, row_sums, out=np.zeros_like(W), where=row_sums != 0)

    loss = np.zeros(len(nodes))

    if source not in idx:
        raise ValueError(f"Source node {source} not found in network nodes: {nodes}")

    source_vec = np.zeros(len(nodes))
    source_vec[idx[source]] = alpha

    for _ in range(n_steps):
        propagated = beta * W_norm.dot(loss)
        loss = np.maximum(source_vec, propagated)
        loss = np.clip(loss, 0, 1)

    return pd.DataFrame({
        "node": nodes,
        "loss": loss,
        "alpha": alpha,
        "beta": beta,
        "n_steps": n_steps,
    })


def estimate_relative_fitness(vulnerability_score, fam50a_loss):
    """
    Convert vulnerability and FAM50A perturbation into a bounded relative fitness score.

    fitness = exp(-vulnerability * FAM50A_loss)

    This makes fitness lower when predicted vulnerability is higher
    and knockdown intensity is stronger.
    """
    return np.exp(-vulnerability_score * fam50a_loss)


alphas = np.linspace(0, 1, 21)
beta = 0.30
n_steps = 5

simulation_rows = []
node_loss_rows = []

for alpha in alphas:
    node_loss = propagate_knockdown(
        G,
        source="FAM50A",
        alpha=alpha,
        beta=beta,
        n_steps=n_steps,
    )

    node_loss_rows.append(node_loss)

    fam50a_loss = float(
        node_loss.loc[node_loss["node"].eq("FAM50A"), "loss"].iloc[0]
    )

    tmp = case_df.copy()
    tmp["alpha"] = alpha
    tmp["fam50a_network_loss"] = fam50a_loss
    tmp["relative_fitness"] = estimate_relative_fitness(
        tmp["vulnerability_positive"],
        fam50a_loss,
    )

    simulation_rows.append(tmp)

simulation = pd.concat(simulation_rows, ignore_index=True)
node_loss_all = pd.concat(node_loss_rows, ignore_index=True)

# ---------------------------------------------------------------------
# Summarize by vulnerability group
# ---------------------------------------------------------------------

summary = (
    simulation
    .groupby(["alpha", "FAM50A_predicted_vulnerability_group"])
    .agg(
        n=("case_submitter_id", "nunique"),
        median_relative_fitness=("relative_fitness", "median"),
        mean_relative_fitness=("relative_fitness", "mean"),
        q25_relative_fitness=("relative_fitness", lambda x: x.quantile(0.25)),
        q75_relative_fitness=("relative_fitness", lambda x: x.quantile(0.75)),
        median_vulnerability=("vulnerability_positive", "median"),
    )
    .reset_index()
)

out_sim = TABLE_DIR / "FAM50A_knockdown_simulation_case_level.csv"
out_summary = TABLE_DIR / "FAM50A_knockdown_simulation_group_summary.csv"
out_node_loss = TABLE_DIR / "FAM50A_knockdown_network_node_loss.csv"

simulation.to_csv(out_sim, index=False, encoding="utf-8-sig")
summary.to_csv(out_summary, index=False, encoding="utf-8-sig")
node_loss_all.to_csv(out_node_loss, index=False, encoding="utf-8-sig")

print(f"\nSaved case-level simulation to: {out_sim}")
print(f"Saved group summary to: {out_summary}")
print(f"Saved node loss table to: {out_node_loss}")

print("\nSimulation summary at alpha = 1:")
print(
    summary[summary["alpha"].eq(1.0)][
        [
            "FAM50A_predicted_vulnerability_group",
            "n",
            "median_vulnerability",
            "median_relative_fitness",
            "q25_relative_fitness",
            "q75_relative_fitness",
        ]
    ].to_string(index=False)
)

# ---------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------

group_order = [
    "high_predicted_vulnerability",
    "intermediate",
    "low_predicted_vulnerability",
]

# Figure 1: fitness trajectory by group
plt.figure(figsize=(7, 5))

for group in group_order:
    sub = summary[summary["FAM50A_predicted_vulnerability_group"].eq(group)]
    if sub.empty:
        continue

    plt.plot(
        sub["alpha"],
        sub["median_relative_fitness"],
        marker="o",
        label=group,
    )

plt.xlabel("FAM50A knockdown intensity alpha")
plt.ylabel("Median relative fitness")
plt.title("Simulated FAM50A knockdown response by predicted vulnerability group")
plt.legend(fontsize=8)
plt.tight_layout()

fig1 = FIG_DIR / "FAM50A_knockdown_simulation_fitness_trajectory.png"
plt.savefig(fig1, dpi=300)

# Figure 2: node loss trajectory
plt.figure(figsize=(7, 5))

for node in node_loss_all["node"].unique():
    sub = node_loss_all[node_loss_all["node"].eq(node)]
    plt.plot(
        sub["alpha"],
        sub["loss"],
        marker="o",
        label=node,
    )

plt.xlabel("FAM50A knockdown intensity alpha")
plt.ylabel("Network perturbation loss")
plt.title("Network propagation of FAM50A knockdown")
plt.legend(fontsize=8)
plt.tight_layout()

fig2 = FIG_DIR / "FAM50A_knockdown_network_node_loss.png"
plt.savefig(fig2, dpi=300)

# Figure 3: relative fitness distribution at full knockdown
full = simulation[simulation["alpha"].eq(1.0)].copy()

data_by_group = [
    full[full["FAM50A_predicted_vulnerability_group"].eq(group)]["relative_fitness"]
    for group in group_order
]

plt.figure(figsize=(7, 5))
plt.boxplot(
    data_by_group,
    tick_labels=[
        "High predicted\nvulnerability",
        "Intermediate",
        "Low predicted\nvulnerability",
    ],
    showfliers=True,
)
plt.ylabel("Relative fitness at alpha = 1")
plt.title("Simulated full FAM50A knockdown effect")
plt.tight_layout()

fig3 = FIG_DIR / "FAM50A_knockdown_alpha1_fitness_boxplot.png"
plt.savefig(fig3, dpi=300)

print("\nSaved figures:")
print(f"  {fig1}")
print(f"  {fig2}")
print(f"  {fig3}")

print("\nDone.")