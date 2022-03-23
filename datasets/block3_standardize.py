# the following fields are meant to serve as a template
control = ???
replace_dict = {
    control: 'control',
}
adata.obs['perturbation_name'] = adata.obs[???].replace(replace_dict)
# remember to use + for combinations
adata.obs['perturbation_type'] = # 'small molecule' or 'genetic'
adata.obs['perturbation_value'] = adata.obs[???]
adata.obs['perturbation_unit'] = # 'ug', 'mg', 'hrs', etc.
