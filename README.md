# CRC Dependency Analysis: FAM50A and FAM50B

This repository contains an early-stage computational analysis of colorectal cancer gene dependency using DepMap CRISPR gene effect data and multi-omics features.

The main result so far is that FAM50A dependency in colorectal cancer cell lines is strongly associated with FAM50B expression.

## Project goal

The broader goal of this project is to build a pipeline for predicting cancer cell gene dependency and later connect the prediction results to CRISPR knockdown network simulation.

Current stage:

- DepMap colorectal cancer cell line analysis
- CRISPR gene dependency prediction
- FAM50A and FAM50B relationship validation

Future stage:

- STRING and Reactome network analysis
- TCGA COAD and READ patient-level extension
- CRISPR knockdown propagation simulation

## Dataset

Data source: DepMap

Cancer type: Colorectal adenocarcinoma

Selected cell lines: 59

Input features:

- Gene expression
- Copy number variation
- Somatic mutation

Prediction label:

- CRISPR gene effect score

Raw DepMap files are not included in this repository because of file size and data-use considerations.

## Analysis workflow

The analysis was performed in the following order.

1. DepMap files were downloaded manually.
2. Colorectal cancer cell lines were selected using Model metadata.
3. CRISPR gene effect, expression, CNV, and mutation data were filtered to 59 colorectal cancer cell lines.
4. Mutation data were converted into a gene-level binary mutation matrix.
5. Initial colorectal cancer genes such as KRAS, BRAF, APC, TP53, PIK3CA, SMAD4, and FBXW7 were tested as dependency prediction targets.
6. Candidate dependency targets were selected based on CRISPR gene effect variability.
7. Elastic Net models were trained using multi-omics features.
8. FAM50A was selected as the strongest pilot dependency target.
9. FAM50B expression was identified as the key explanatory feature for FAM50A dependency.

## Key result 1: FAM50A dependency prediction

FAM50A was the strongest pilot target among the tested candidate dependency genes.

Full multi-omics model performance:

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
- Full model: R2 = 0.404
- Full model without FAM50B expression: R2 = 0.103

This suggests that FAM50B expression is the main explanatory feature for FAM50A dependency in the selected colorectal cancer cell lines.

## Key result 3: FAM50B-low cells show stronger FAM50A dependency

FAM50B expression was split into low and high groups by median expression.

Group comparison:

- FAM50B-low group median FAM50A gene effect: -1.36
- FAM50B-high group median FAM50A gene effect: -0.49
- Mann-Whitney p-value: 1.06e-05

Because more negative CRISPR gene effect values indicate stronger dependency, this result suggests that FAM50B-low colorectal cancer cell lines may be more dependent on FAM50A.

## Interpretation

The current result suggests that FAM50A may act as a context-specific vulnerability in FAM50B-low colorectal cancer cell lines.

This should be interpreted as a computational association rather than direct experimental proof.

## Output files

Important result tables are stored in artifacts/tables.

Important figures are stored in artifacts/figures.

Representative outputs include:

- FAM50A prediction diagnostics
- FAM50A permutation test
- FAM50B expression versus FAM50A dependency plot
- FAM50A dependency comparison between FAM50B-low and FAM50B-high groups

## Next steps

1. Build a STRING functional network using FAM50A, FAM50B, and stable model features.
2. Perform Reactome pathway enrichment analysis.
3. Extend the analysis to TCGA COAD and READ patient tumors.
4. Define FAM50B-low patient groups.
5. Develop a CRISPR knockdown network propagation simulation.
