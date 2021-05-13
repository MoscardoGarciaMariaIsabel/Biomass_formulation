Folder content:

tsv files containing data matrices (described below)
for 17,995  genes x 325 samples (= 324 cell lines + HT29v1.1).
 
- 00_logFCs.tsv: depletion log fold changes
- 01_corrected_logFCs.tsv: CRISPRcleaned depletion log fold changes
- 01a_qnorm_corrected_logFCs.tsv: quantile normalized CRISPRcleaned depletion log fold changes
- 02a_BayesianFactors.tsv: Bayesian factors (obtained with BAGELr from 01).
- 02b_MageckFDRs.tsv: Mageck depletion FDRs (obtained with Mageck from corrected sgRNA counts derived from 01).
- 03_scaledBayesianFactors.tsv: Scaled Bayesian Factors (obtained from 02a)
- 04_binaryDepScores.tsv: Binary depletion scores (derived from 03).