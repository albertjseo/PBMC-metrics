import pandas as pd
import scanpy as sc

# Set default matplotlib params
sc.settings.verbosity = 0  # providing additional information / feedback from program. 0 - errors, 1 - warning messages, 2 - informational message, 3 - hints, 4 - debugging
sc.logging.print_header() #Versions that might influence the numerical results. Matplotlib and Seaborn are excluded from this.
sc.settings.set_figure_params(dpi=100, facecolor="white") # sets resolutions, sizing, styling, and formatting of figure.

results_file = "write/pbmc3k.h5ad"  # the file that will store the analysis results

# adata is an AnnData object that can be sliced like a dataframe. AnnData stores a data matrix with observations, variables, and unstructures annotations.
# read_10x_mtx returns an AnnData object
adata = sc.read_10x_mtx(
    "/Users/albertseo/Development/PBMC-metrics/filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names
    cache=True,  # write a cache file for faster subsequent reading
)

adata.var_names_make_unique() # apply common gene name

# Preprocessing
sc.pl.highest_expr_genes(adata, n_top=50) # Plotting top 50 expressing genes from adata
sc.pp.filter_cells(adata, min_genes=200) # Filter cell outliers based on counts and numbers of genes expressed.
sc.pp.filter_genes(adata, min_cells=3) # Filter genes based on number of cells or counts.

# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True) # Returns calculated metrics as a dataframe.

# Slice the  to filter it
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data to complete the normalization
sc.pp.log1p(adata)

# Identify high variance genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Plot high variance genes (darker) transposed ontop of non-high variance genes (lighter)
sc.pl.highly_variable_genes(adata)

# Principal Component Analysis
# tl is a tool that is used to transform data that is not preprocessing
sc.tl.pca(adata, svd_solver="arpack") # PCA applies principal component analysis

# plot the PCA and color by gene expression (louvain would result in Louvain clusters)
sc.pl.pca(adata, color="CST3")
