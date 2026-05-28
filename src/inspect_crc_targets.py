from pathlib import Path
import pandas as pd

INTERIM_DIR = Path("data/interim")

gene_effect = pd.read_csv(INTERIM_DIR / "depmap_crc_gene_effect.csv", nrows=0)
expression = pd.read_csv(INTERIM_DIR / "depmap_crc_expression.csv", nrows=0)
cnv = pd.read_csv(INTERIM_DIR / "depmap_crc_cnv.csv", nrows=0)
mutation = pd.read_csv(INTERIM_DIR / "depmap_crc_mutation_matrix.csv", nrows=0)

targets = ["KRAS", "BRAF", "PIK3CA", "APC", "TP53", "SMAD4", "FBXW7"]

def find_gene_columns(columns, gene):
    hits = []
    for col in columns:
        col_str = str(col)
        if col_str == gene or col_str.startswith(gene + " ") or col_str.startswith(gene + " ("):
            hits.append(col_str)
    return hits

print("Checking target genes...\n")

for gene in targets:
    print(f"=== {gene} ===")

    ge_hits = find_gene_columns(gene_effect.columns, gene)
    expr_hits = find_gene_columns(expression.columns, gene)
    cnv_hits = find_gene_columns(cnv.columns, gene)
    mut_col = f"mut_{gene}"

    print(f"gene_effect: {ge_hits[:5]}")
    print(f"expression:  {expr_hits[:5]}")
    print(f"cnv:         {cnv_hits[:5]}")
    print(f"mutation:    {'YES' if mut_col in mutation.columns else 'NO'}")
    print()