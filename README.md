# CRC Dependency Analysis: FAM50A and FAM50B

This repository contains a computational pipeline for predicting colorectal cancer gene dependency and connecting the result to TCGA patient stratification and CRISPR knockdown network simulation.

The main finding so far is that FAM50A dependency in colorectal cancer cell lines is strongly associated with FAM50B expression. This relationship was detected in DepMap cell lines, projected onto TCGA COAD/READ tumors, and used to build a simple knockdown simulation model.

## Project goal

The broader goal of this project is to build an AI-based workflow for cancer cell gene dependency prediction and CRISPR knockdown network simulation.

Current completed modules:

- DepMap colorectal cancer cohort construction
- Multi-omics dependency prediction
- FAM50A and FAM50B relationship validation
- STRING and Reactome network/pathway analysis
- TCGA COAD/READ patient projection
- Clinical association analysis
- CRISPR knockdown network simulation MVP

## Data sources

The analysis uses public cancer genomics and functional genomics resources.

DepMap cell line data:

- CRISPR gene effect scores
- Gene expression
- Copy number variation
- Somatic mutation
- Model metadata

TCGA/GDC patient data:

- TCGA-COAD and TCGA-READ STAR-Counts expression files
- Clinical metadata including stage, vital status, survival follow-up, age, and project information

Raw DepMap and GDC files are not included in this repository because of file size and data-use considerations. The repository stores analysis scripts, summary tables, figures, and manifests.

## Analysis workflow

The project was developed in the following order.

1. DepMap files were downloaded and checked.
2. Colorectal adenocarcinoma cell lines were selected from Model metadata.
3. CRISPR gene effect, expression, CNV, and mutation matrices were filtered to the colorectal cancer cohort.
4. Mutation data were converted into a gene-level binary mutation matrix.
5. Initial colorectal cancer genes such as KRAS, BRAF, APC, TP53, PIK3CA, SMAD4, and FBXW7 were tested as dependency targets.
6. Candidate dependency targets were selected based on CRISPR gene effect variability.
7. Elastic Net models were trained using expression, CNV, and mutation features.
8. FAM50A was selected as the strongest pilot dependency target.
9. FAM50B expression was identified as the key explanatory feature for FAM50A dependency.
10. STRING and Reactome analyses were performed for model-derived genes.
11. TCGA COAD/READ expression data were used to project FAM50A vulnerability into patient tumors.
12. Clinical association analyses were performed using stage, survival, age, and project metadata.
13. A simple FAM50A knockdown network simulation was implemented as an MVP model.

## Key result 1: FAM50A dependency prediction in DepMap

FAM50A was the strongest pilot target among candidate dependency genes.

Full multi-omics Elastic Net model:

- R2: 0.404
- Correlation: 0.648
- RMSE: 0.488

Permutation test:

- Observed R2: 0.404
- Maximum permuted R2: 0.066
- Empirical p-value: 0.0196

This suggests that the FAM50A dependency prediction signal is unlikely to be explained by random label structure alone.

## Key result 2: FAM50B expression explains FAM50A dependency

FAM50B expression was the most stable feature selected across cross-validation folds.

Model comparison:

- FAM50B expression only: R2 = 0.434
- Full multi-omics model: R2 = 0.404
- Full model without FAM50B expression: R2 = 0.103

This suggests that FAM50B expression is the main explanatory variable for FAM50A dependency in the selected colorectal cancer cell lines.

## Key result 3: FAM50B-low cell lines show stronger FAM50A dependency

FAM50B expression was split into low and high groups by median expression.

Group comparison:

- FAM50B-low group median FAM50A gene effect: -1.36
- FAM50B-high group median FAM50A gene effect: -0.49
- Mann-Whitney p-value: 1.06e-05

Because more negative CRISPR gene effect values indicate stronger dependency, this result suggests that FAM50B-low colorectal cancer cell lines may be more dependent on FAM50A.

## Network and pathway analysis

STRING functional network analysis returned a compact network connecting FAM50A, FAM50B, and ZDBF2.

Network structure:

- FAM50A -- FAM50B
- FAM50B -- ZDBF2

This supports the interpretation that FAM50A and FAM50B are functionally associated, although STRING evidence should not be treated as direct experimental proof of physical interaction.

Reactome over-representation analysis was also performed using model-derived gene sets. No pathway passed the conventional FDR < 0.05 threshold. Weak exploratory signals appeared for intestinal absorption and mitochondrial RNA modification, but these should not be interpreted as statistically significant pathway-level findings.

## TCGA patient projection

TCGA COAD/READ primary tumor expression data were used to project the DepMap-derived FAM50A vulnerability relationship into patient tumors.

Expression cohort:

- TCGA primary tumor samples: 647
- TCGA-COAD: 481
- TCGA-READ: 166

FAM50B expression groups:

- FAM50B-low tumors: 324
- FAM50B-high tumors: 323
- FAM50B median log2(TPM+1): 2.4094

FAM50A expression was slightly higher in the FAM50B-high group than in the FAM50B-low group.

FAM50A and FAM50B expression correlation in TCGA:

- Pearson r: 0.224
- Spearman r: 0.174

This indicates that TCGA tumors can be stratified by FAM50B expression, but TCGA expression alone does not directly measure CRISPR dependency.

## Predicted FAM50A vulnerability in TCGA

A simple DepMap-derived linear model was trained using FAM50B expression to predict FAM50A gene effect.

DepMap FAM50B-only model:

- Intercept: -1.3587
- Slope: 0.3651
- R2: 0.5179
- RMSE: 0.4385

The model was applied to TCGA COAD/READ tumors to generate a predicted FAM50A vulnerability score.

TCGA predicted vulnerability groups:

- High predicted vulnerability: 162 tumors
- Intermediate: 323 tumors
- Low predicted vulnerability: 162 tumors

Project-level comparison:

- COAD median vulnerability score: 0.5591
- READ median vulnerability score: 0.2303
- Mann-Whitney p-value: 0.001198

This suggests that predicted FAM50A vulnerability is more common or stronger in TCGA-COAD than TCGA-READ in this model.

## Clinical association analysis

Predicted FAM50A vulnerability scores were merged with TCGA clinical metadata.

Merged cohort:

- Merged cases: 624
- Stage-known cases: 506
- Survival-usable cases: 530

Stage association:

- Stage-vulnerability chi-square p-value: 0.2733

Survival analysis:

- High vs low vulnerability log-rank p-value: 0.09184
- High vs low overall survival days Mann-Whitney p-value: 0.2737

Cox regression:

- Vulnerability-only HR: 0.933, p = 0.436
- Project-adjusted HR: 0.928, p = 0.402
- Project and age-adjusted HR: 0.891, p = 0.205
- Project, stage, and age-adjusted HR: 0.961, p = 0.721

These results do not support predicted FAM50A vulnerability as an independent overall survival prognostic marker in the current TCGA analysis. The safer interpretation is that FAM50A vulnerability is a molecular dependency hypothesis rather than a validated survival marker.

## CRISPR knockdown network simulation MVP

A simple FAM50A knockdown simulation was implemented using the TCGA predicted vulnerability score and the STRING-derived network.

Simulation setup:

- Source node: FAM50A
- Network: FAM50A -- FAM50B -- ZDBF2
- Knockdown intensity alpha: 0 to 1
- Relative fitness model: fitness = exp(-vulnerability score x FAM50A loss)

Simulated full knockdown result at alpha = 1:

- High predicted vulnerability group median relative fitness: 0.385
- Intermediate group median relative fitness: 0.619
- Low predicted vulnerability group median relative fitness: 1.000

This shows that the predicted high-vulnerability group has the strongest simulated fitness decrease under FAM50A knockdown.

This simulation is a conceptual MVP model. It does not replace wet-lab CRISPR knockdown experiments or mechanistic biochemical modeling.

## Main interpretation

The current pipeline supports the following interpretation.

FAM50B-low colorectal cancer cell lines show stronger FAM50A dependency in DepMap. This relationship can be projected onto TCGA COAD/READ tumors to define a predicted FAM50A vulnerability axis. The axis does not appear to be an independent overall survival marker in the current TCGA clinical analysis, but it can be used to define a molecularly stratified group that may be more sensitive to FAM50A perturbation in a simulation setting.

## Limitations

Important limitations:

- The DepMap cohort contains only 59 colorectal cancer cell lines.
- TCGA tumors do not have direct CRISPR dependency measurements.
- The TCGA predicted vulnerability score is inferred from a DepMap-derived expression model.
- Clinical association results are exploratory and affected by missing clinical metadata.
- Reactome pathway enrichment did not identify statistically significant pathways at FDR < 0.05.
- The knockdown network simulation is a toy MVP model, not a validated mechanistic model.
- Experimental validation would be required to confirm FAM50A as a therapeutic vulnerability in FAM50B-low colorectal cancer.

## Repository structure

Main folders:

- src: analysis scripts
- artifacts/tables: result tables and summaries
- artifacts/figures: generated plots
- artifacts/manifests: file manifests and metadata tables
- data: local raw, interim, and processed data folders excluded from GitHub

## Representative outputs

Representative figures include:

- FAM50A actual vs predicted dependency plot
- FAM50A permutation test
- FAM50B expression versus FAM50A dependency plot
- FAM50A dependency by FAM50B group
- STRING network around FAM50A model features
- TCGA FAM50B expression distribution
- TCGA predicted FAM50A vulnerability distribution
- Kaplan-Meier curve by predicted vulnerability group
- FAM50A knockdown simulation trajectory

Representative tables include:

- Candidate dependency target results
- FAM50A ablation results
- FAM50B feature importance validation
- Reactome enrichment results
- TCGA FAM50 expression summary
- TCGA predicted vulnerability summary
- TCGA clinical association summary
- Cox regression summary
- Knockdown simulation summary

## Next steps

Possible next steps:

1. Refine the knockdown simulation model using a larger biological network.
2. Add confidence intervals or bootstrapping to the TCGA vulnerability projection.
3. Compare FAM50A/FAM50B patterns across additional cancer types.
4. Add independent cell line or CRISPR screen validation if available.
5. Organize the repository into modular folders for preprocessing, modeling, TCGA analysis, and simulation.
6. Write a concise research report based on the current pipeline and results.
