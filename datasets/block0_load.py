author_year = #e.g. You_2022
is_counts = # True/False
var_genes = # field in .var or None
doi = # of paper source

import numpy as np
import pandas as pd
import scanpy as sc
sc.set_figure_params(dpi=100, frameon=False)
sc.logging.print_header()

# verify
assert(doi in pd.read_csv('../personal.csv').DOI.values)

adata = sc.read(f'{author_year}_raw.h5ad')
adata
