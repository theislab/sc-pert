# metalabels
adata.uns['preprocessing_nb_link'] = f'https://nbviewer.org/github/theislab/sc-pert/blob/main/datasets/{author_year}_curation.ipynb'
adata.uns['doi'] = doi
display(adata.obs.describe(include='all').T.head(20))

# use gene symbols as gene names
if var_genes:
    adata.var[var_genes] = adata.var[var_genes].astype(str)
    adata.var = adata.var.reset_index().set_index(var_genes)
    print(adata.var_names)

# filtering and processing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=20)
adata.var_names_make_unique()

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
