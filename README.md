\# CRC Dependency Analysis: FAM50A and FAM50B



\## Project Summary



This project analyzes colorectal cancer cell lines from DepMap to predict CRISPR-based gene dependency using multi-omics features.



The current result focuses on a candidate dependency relationship between FAM50A and FAM50B expression.



\## Dataset



| Item | Description |

|---|---|

| Data source | DepMap |

| Cancer type | Colorectal adenocarcinoma |

| Number of cell lines | 59 |

| Features | Expression, CNV, mutation |

| Label | CRISPR gene effect score |



Raw DepMap files are not included in this repository.



\## Workflow



1\. DepMap files were downloaded manually.

2\. Colorectal cancer cell lines were selected using Model metadata.

3\. CRISPR gene effect, expression, CNV, and mutation data were filtered to 59 colorectal cancer cell lines.

4\. Mutation data were converted into a gene-level binary matrix.

5\. Initial colorectal cancer genes such as KRAS, BRAF, APC, TP53, PIK3CA, SMAD4, and FBXW7 were tested.

6\. Candidate dependency targets were selected based on CRISPR gene effect variability.

7\. Elastic Net models were trained using multi-omics features.

8\. FAM50A was selected as the strongest pilot dependency target.

9\. FAM50B expression was identified as the key explanatory feature.



\## Key Results



\### FAM50A Dependency Prediction



| Model | Result |

|---|---:|

| Full multi-omics model R² | 0.404 |

| Correlation | 0.648 |

| Permutation test p-value | 0.0196 |



The observed prediction performance was higher than all 50 permuted-label models.



\## FAM50B Expression as the Main Feature



| Feature set | R² |

|---|---:|

| FAM50B expression only | 0.434 |

| Full model | 0.404 |

| Full model without FAM50B expression | 0.103 |



This suggests that FAM50B expression is the main explanatory feature for FAM50A dependency.



\## Group Comparison



| Group | Median FAM50A gene effect |

|---|---:|

| FAM50B-low group | -1.36 |

| FAM50B-high group | -0.49 |



Mann-Whitney p-value: 1.06e-05



Because more negative CRISPR gene effect values indicate stronger dependency, this result suggests that FAM50B-low colorectal cancer cell lines may be more dependent on FAM50A.



\## Interpretation



The current result suggests that FAM50A may act as a context-specific vulnerability in FAM50B-low colorectal cancer cell lines.



This is a computational association and should not be interpreted as direct experimental proof.



\## Next Steps



1\. Build a STRING network using FAM50A, FAM50B, and stable model features.

2\. Perform Reactome pathway enrichment analysis.

3\. Extend the analysis to TCGA COAD and READ patient tumors.

4\. Define FAM50B-low patient groups.

5\. Develop a CRISPR knockdown network simulation.

