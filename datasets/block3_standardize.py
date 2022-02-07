adata = sc.read(f'{author_year}.h5ad')
# the following fields are meant to serve as a template
adata.obs['perturbation_name'] = ???
adata.obs['perturbation_type'] = # 'small molecule' or 'genetic'
adata.obs['perturbation_value'] = adata.obs[???]
adata.obs['perturbation_unit'] = # 'ug', 'mg', 'hrs', etc.
