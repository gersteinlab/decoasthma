# decoasthma <img src="man/figures/logo.png" align="right" height="160" alt="Hex sticker for decoasthma: a pair of lungs, a budding yeast, and a bacterial rod." />

[![DOI](https://zenodo.org/badge/255353495.svg)](https://zenodo.org/badge/latestdoi/255353495)

Code to reproduce the analyses in Spakowicz, Lou, et al., *Genome Biology* (2020). The paper introduces LDA-link, a method for relating microbes to genes in heterogeneous or noisy RNA-seq data from the sputum of asthmatic patients.

## Citation

Spakowicz D\*, Lou S\*, Barron B, Gomez JL, Li T, Liu Q, Grant N, Yan X, Hoyd R, Weinstock G, Chupp GL, Gerstein M (2020) Approaches for integrating heterogeneous RNA-seq data reveal cross-talk between microbes and genes in asthmatic patients. *Genome Biology* 21:150. https://doi.org/10.1186/s13059-020-02033-z · [PubMed: 32571363](https://pubmed.ncbi.nlm.nih.gov/32571363/) (\* contributed equally)

A preprint is posted on bioRxiv: https://doi.org/10.1101/765297

The data are being submitted to dbGaP under BioProject SUB7102729 and are available through academic collaboration. Please contact daniel.spakowicz@osumc.edu and shaoke.lou@yale.edu with any questions.

<br clear="all">

## Graphical abstract

<img src="man/figures/graphical-abstract.png" width="100%" alt="Graphical abstract. Induced sputum RNA-seq from 115 asthmatic patients is split into a human gene table and a microbe table. Gene expression is deconvolved into neutrophil, eosinophil, and mast-cell fractions and checked against single-cell RNA-seq and microscopy. Latent Dirichlet allocation summarizes genes and microbes as ten topics each. LDA-link trains a random forest on those topics and reports 1,883 high-confidence links, including Haemophilus with IL1B in mast cells and Candida with GCSAML in eosinophils.">

## Figure scripts

Plotting code for the manuscript figures lives in directories named for each figure. The scripts read processed count matrices, clinical tables, and intermediate `.Rdata` files prepared by the notebooks in [`data-processing/`](data-processing). Those input files are shared through the dbGaP submission above and by academic collaboration.

| Manuscript figure | Script | What it draws |
| --- | --- | --- |
| Figure 1, biotype stacked bars | [`Figure1/Fig1_biotype_stacked-bar.Rmd`](Figure1/Fig1_biotype_stacked-bar.Rmd) | Read-alignment summary by biotype (`fig1_biotype_stacked-bar.pdf`), plus a supplemental stacked bar |
| Figure 2B–C | [`Figure2/B&C/novel-cell-types.Rmd`](Figure2/B%26C/novel-cell-types.Rmd) | Single-cell RNA-seq clusters, reference-cell labels, and per-patient cell-type fractions |
| Figure 2D and 2F | [`Figure2/F/deconvolution-figure.Rmd`](Figure2/F/deconvolution-figure.Rmd) | Cell-fraction estimates from the gene table, comparison with microscopy (cytospin) counts, and correlations with clinical variables. Also writes the NMF comparison used in the supplement |
| Figure 2E | [`Figure2/Fig2.Rmd`](Figure2/Fig2.Rmd) | Correlation of cell-type fractions with the ten LDA gene topics (`lm22-lda10.png`). This notebook also repeats the Figure 2B–D and 2F code |
| Figure 3A–C | [`Figure3/A_B_C/exogenous-figure.Rmd`](Figure3/A_B_C/exogenous-figure.Rmd) | Exogenous (microbial) abundances, associations with clinical variables and cell fractions, and phylum-level clustering by asthma severity |
| Figure 3D and Figure S5 | [`Figure3/Fig3D.r`](Figure3/Fig3D.r) | Microbe co-abundance networks, including the network with an LDA overlay |
| Figure 4A–F and Figures S6–S7 | [`Figure 4-LDA-Link/fig4.LDA_link.r`](Figure%204-LDA-Link/fig4.LDA_link.r) | LDA-link training pairs, pathway enrichment, random-forest topic importance, and topic membership. The Figure 4A density plot sits in an `if (FALSE)` block near the top of the script |
| Figure 5A–B | [`Figure 5/fig5.r`](Figure%205/fig5.r) | Bipartite and tripartite graphs of gene–microbe links, and cell-type expression for linked genes such as IL1B, GCSAML, and CACNA1E |

`Figure2/Fig2.Rmd` is the combined Figure 2 notebook. For a single panel, the shorter files under `Figure2/B&C/` and `Figure2/F/` are the more direct route, except for panel E, which is only in the combined notebook.

Upstream of those figures, [`data-processing/`](data-processing) builds the matrices the plotting scripts expect:

- [`01_clinical-data-subset.Rmd`](data-processing/01_clinical-data-subset.Rmd) subsets the clinical table
- [`02_quantile-normalize.R`](data-processing/02_quantile-normalize.R) quantile-normalizes expression
- [`03_batch-effects-check.Rmd`](data-processing/03_batch-effects-check.Rmd) checks batch effects
- [`generate_downloadKingdom.R`](data-processing/generate_downloadKingdom.R) and [`generate_exceRpt_lsf_sub.sh`](data-processing/generate_exceRpt_lsf_sub.sh) prepare the exceRpt alignments that separate human and exogenous reads

## Running LDA-link

LDA-link is the script [`Figure 4-LDA-Link/fig4.LDA_link.r`](Figure%204-LDA-Link/fig4.LDA_link.r). It fits a 10-topic latent Dirichlet allocation model to the filtered gene counts and to the microbe counts, trains a random forest on the topic weights of gene–microbe pairs, and then scores every pair. Source it from R after the packages below are attached and `wd` at the top of the script points at your data directory. The script calls `setwd(wd)` itself.

```r
install.packages(c("randomForest", "PRROC", "glmnet", "e1071"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("topicmodels")

library(randomForest)
library(PRROC)
library(glmnet)
library(e1071)

# Edit wd inside the script, then from the repository root:
source("Figure 4-LDA-Link/fig4.LDA_link.r")
```

`randomForest()` is called in the cross-validation loop before the script reaches `library(randomForest)`, so attach that package first or the run stops there. `topicmodels` is loaded in the script, just before the gene model is fit. The main path uses random forests (`method = "rf"`). `glmnet` and `e1071` are loaded for the other classifiers in that loop. The function `ctool()` at the bottom of the file is a cross-validation helper and is not called on a normal run; calling it also needs `ROCR` and `AUC`.

### Inputs

`wd` is the directory that contains the expression matrix. The two exogenous files are read from `../exogenous/` relative to `wd`. Create `0rpmlogF_noCtrl/` inside `wd` before sourcing; the importance plot and the link tables are written there. None of these inputs are in the git repository.

| Path relative to `wd` | What the script expects |
| --- | --- |
| `counts.rpm.protein.rpkm.clinical.Rdata` | `all.mats.protein$rpm`, genes by samples, in reads per million |
| `../exogenous/exo.genus.wRowName.rdata` | `exo.genus`, genus-level abundances, samples in the columns |
| `../exogenous/ExoAsthma_humanGenes_clinical.RData` | `df`, the clinical table. Column 1 is the sample id and column 46 is asthma severity |
| `exo_signal2_ldad10.txt` | Tab-separated counts with a header and row names. Rows are samples and columns are microbes. This matrix is the one passed to the microbe LDA |
| `fig4a_corlinks.gene.david.GAD.txt` | DAVID Genetic Association Database table, used only to draw Figure 4B |
| `0rpmlogF_noCtrl/` | Output directory. Create it first |

Sample ids have to line up. Expression colnames have a trailing `.fq` removed, and the script then matches them to the columns of `exo.genus`. Clinical ids in `df[, 1]` are matched to LDA document names after `.` in those names is replaced with `-`. The columns of `exo_signal2_ldad10.txt` have to be the microbes kept in the filtered genus table (present in more than 10 samples), in that same row order: pair indexes are shared between the correlation matrix and the LDA topic matrix.

### What a run does

1. Drops microbes with no reads, replaces missing abundances with zero, and keeps genes that are above the per-sample median in more than 30 samples (`cutoff = 30`). A stricter count, `bulk.filt` at 70 samples, is computed and is not the matrix passed to LDA.
2. Fits a 10-topic Gibbs LDA to the genes with `topicmodels::LDA`. Counts are `ceiling(RPM / 10)`, capped at 1000. Saves `bulk_ldaout10.rdata`.
3. Fits a 10-topic Gibbs LDA to `exo_signal2_ldad10.txt` with seed 123, burn-in 1000, thinning 100, and 1000 iterations. Saves `exo_ldaout10.rdata`.
4. Writes `top20_bulk_gene2topic_dist.pdf` and `top10_exo_microbe2topic_dist.pdf` (Figure 4D–F and Figures S6–S7).
5. Builds the training set. Positive pairs have absolute Pearson correlation above `cor.cut` (0.4). Negative pairs have absolute correlation below 0.05 and are sampled down to the same number of pairs. The feature vector for each pair is the 10 gene-topic weights followed by the 10 microbe-topic weights. Adjusted p-values are stored in `be.cor.padj_fdr`. The published labels also required p < 1e−5 and FDR < 0.016; this script applies the correlation cutoffs only.
6. Runs 10-fold cross-validation with a 500-tree random forest and prints the ROC AUC and the precision-recall AUC (`roc1`, `prc1`).
7. Trains a final forest, `rf.b2e`, on all labeled pairs and writes `0rpmlogF_noCtrl/rf_varimp.pdf`.
8. Scores every gene–microbe pair. The prediction is split at 1.2 million rows. Pairs with probability above `prob.cut` are written out.

The manuscript calls a link at probability > 0.95. In the script `prob.cut` is `0.9`, with `0.95` noted in a comment beside it. Set `prob.cut <- 0.95` before the `write.table` calls to use the published cutoff. The link tables are:

- `0rpmlogF_noCtrl/b2e.links.probcut<cut>_corcut0.8.rpm.las1se.txt`, gene and microbe, no header
- `0rpmlogF_noCtrl/b2e.links.c<cut>_corcut0.8.rpm.las1se.gene.txt`, the unique genes
- `rf.link.2ndRun.rdata`, the full workspace, including the probability matrix `b2e.prob`

Figure 4A, the correlation-density plot, sits in `if (FALSE)` near the top of the script. Change that condition to `TRUE` to write `bulk_exo_cor_distr.pdf`. Figure 4B is `fig4B.pdf`, drawn from the DAVID table. Figure 4C is `Fig4c.pdf`, from `varImpPlot(rfmodel)`, which is the forest fit on the last cross-validation fold. The importance plot for the final model is `rf_varimp.pdf`.

[`Figure 5/fig5.r`](Figure%205/fig5.r) reads its own link table, `b2e.links_rf.txt` (no header; column 1 is the gene and column 2 is the microbe), plus `fig5_pin2gene.rdata`, both from the `Figure 5` directory.

The microbe LDA is seeded. The gene LDA and the sampling of negative pairs and cross-validation folds are not, so the called links can differ from run to run. Scoring every pair is the long step.
