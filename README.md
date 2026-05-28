\# AI-based Prediction of Cancer Cell Gene Dependency and CRISPR Knockdown Network Simulation



\## Overview



This project explores gene dependency prediction in colorectal cancer cell lines using DepMap CRISPR gene effect scores and multi-omics features. The current analysis focuses on DepMap colorectal adenocarcinoma cell lines and identifies a candidate context-specific vulnerability involving FAM50A dependency and FAM50B expression.



\## Current Dataset



\* Data source: DepMap

\* Cancer type: Colorectal adenocarcinoma

\* Number of selected cell lines: 59

\* Features used:



&#x20; \* Gene expression

&#x20; \* Copy number variation

&#x20; \* Somatic mutation

\* Label:



&#x20; \* CRISPR gene effect score



Raw DepMap data files are not included in this repository. Only analysis scripts, summary tables, and figures are stored.



\## Workflow Completed So Far



1\. DepMap files were manually downloaded.

2\. Colorectal cancer cell lines were filtered using Model metadata.

3\. CRISPR gene effect, expression, CNV, and mutation matrices were subset to the colorectal cancer cohort.

4\. Mutation data were transformed into a gene-level binary mutation matrix.

5\. Initial colorectal cancer biomarker genes such as KRAS, BRAF, APC, TP53, PIK3CA, SMAD4, and FBXW7 were tested as dependency prediction targets.

6\. Candidate target genes were selected based on CRISPR gene effect variability.

7\. Elastic Net models were trained using expression, CNV, and mutation features.

8\. FAM50A was selected as the strongest pilot dependency target.

9\. FAM50B expression was identified as the key explanatory feature for FAM50A dependency.



\## Key Results



\### FAM50A Dependency Prediction



The best pilot target was FAM50A.



\* Full multi-omics model:



&#x20; \* R² = 0.404

&#x20; \* Correlation = 0.648



\* Permutation test:



&#x20; \* Observed R² = 0.404

&#x20; \* Maximum permuted R² = 0.066

&#x20; \* Empirical p-value = 0.0196



These results suggest that the FAM50A dependency prediction signal is unlikely to be explained by random label structure alone.



\### FAM50B Expression as a Key Feature



FAM50B expression was the most stable feature selected across cross-validation folds.



Model comparison:



\* FAM50B expression only:



&#x20; \* R² = 0.434

\* Full model:



&#x20; \* R² = 0.404

\* Full model without FAM50B expression:



&#x20; \* R² = 0.103



This suggests that FAM50B expression is the major explanatory variable for FAM50A dependency in the selected colorectal cancer cell lines.



\### Group Comparison



FAM50B expression was split into low and high groups by median expression.



\* FAM50B-low group:



&#x20; \* Median FAM50A gene effect = -1.36

\* FAM50B-high group:



&#x20; \* Median FAM50A gene effect = -0.49

\* Mann-Whitney p-value = 1.06e-05



Since more negative CRISPR gene effect values indicate stronger dependency, this result suggests that FAM50B-low colorectal cancer cell lines may be more dependent on FAM50A.



\## Interpretation



The current result supports the hypothesis that FAM50A may act as a context-specific vulnerability in FAM50B-low colorectal cancer cell lines. This should be interpreted as a computational association rather than direct experimental proof.



\## Next Steps



1\. Build a STRING functional network using FAM50A, FAM50B, and stable model features.

2\. Perform Reactome pathway enrichment analysis.

3\. Extend the analysis to TCGA COAD/READ patient tumors.

4\. Define FAM50B-low patient groups and estimate potential FAM50A vulnerability.

5\. Develop a CRISPR knockdown network propagation simulation.



