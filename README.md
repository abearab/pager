# pager
Container of codes and data to run PAGE algorithm and enable further data curations for hypothesis generation or testing.

## Installation

**Requirements:**

Follw the instructions to install iPAGE and TEISER from the following links:
- https://tavazoielab.c2b2.columbia.edu/iPAGE/ | [GitHub](https://github.com/hanig/PAGE)
- https://tavazoielab.c2b2.columbia.edu/TEISER/ | [GitHub](https://github.com/hanig/TEISER)

**Conda environment:**

Create a conda environment and set the following environment variables to the paths of `pager` and `TEISER` directories.
```bash
conda env create -n pager
conda activate pager
conda env config vars set PAGEDIR="/path/to/PAGE"
conda env config vars set TEISERDIR="/path/to/TEISER"

conda deactivate  # reactivate 
conda activate pager
```

## Background and usage
Codes in this repo generate automated plots and result tables using orignial tools but it enables further data cleaning and imputation to keep / assemble biological meaningful results for hypothesis generation or testing. A python reimplementation of the PAGE algorithm is also provided at https://github.com/goodarzilab/pypage. However, this repo is focused on the original tools mostly implented in Perl.

### PAGE algorithm
[PAGE](https://github.com/hanig/PAGE) algorithm is a set enrichment analysis tool. It's orignially developed by @hanig for analysis of gene expression but it can be generalized for other data types.

Briefly, PAGE quantizes differential measurements into equally populated bins and then, for every given geneset, calculates the mutual information (MI) between each cluster bin and a binary vector of pathway memberships for genes in a given gene set. The significance of each MI value is then assessed through a randomization-based statistical test and hypergeometric distribution to determine whether there is over or under representation of a gene set in each cluster bin. The final result is a p-values matrix in which rows are gene sets and columns are cluster bins (visualized as heatmaps).

### iPAGE run for gene sets
The iPAGE algorithm is a useful tool for gene set and pathway enrichment analysis, e.g. on differential RNA expression values. Here, MSigDB (version 7.4.) gene sets are downloaded and modified to be compatible with iPAGE workflow, see `annotation` directory.

iPAGE in continuous mode accepts gene-level numeric values (e.g., logFCs) as input and performs a PAGE run for each gene set. The output is a p-values matrix where each row corresponds to a gene set and each column corresponds to a cluster bin. The p-values matrix is then used to generate heatmaps for each gene set.

### onePAGE run for single gene set
Based on the biological question, onePAGE can be used to run PAGE algorithm for a single gene set. This is useful when the user is interested in a specific gene set and wants to see how it is distributed across the cluster bins. This analysis can be performed for multiple inputs (e.g., differentially expressed genes from different conditions) and the results can be merged to show the results side by side as heatmaps.

For selected gene sets, the PAGE run is performed using the `run_mi_gene_list.pl` module from TEISERv1.1 with gene-level numeric values (e.g., log fold changes) as input. Then, the p-values matrices of a single gene set for several inputs are merged together to form a single matrix where each row corresponds to an input and each column corresponds to a cluster bin. Finally, heatmaps were generated using the `teiser_draw_matrix.pl` module from TEISERv1.1. The code for TEISERv1.1 can be found at https://github.com/hanig/TEISER/.