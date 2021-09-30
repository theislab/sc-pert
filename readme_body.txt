# sc-pert - Machine learning for perturbational single-cell omics

*This repository provides a community-maintained summary of models and datasets. It was initially curated for [(Cell Systems, 2021)](https://doi.org/10.1016/j.cels.2021.05.016).*

### Objectives

There are various resources for evaluation of single cell perturbation models. We discuss five tasks in the publication which can be supported by the following annotations:

- [GDSC](https://www.cancerrxgene.org/downloads/bulk_download) provides a collection of cell viability measurements for many compounds and cell lines. We provide a [code snippet](https://github.com/theislab/sc-pert/blob/main/resources.py#L4) to conveniently load GDSC-provided z-score compound response rankings per cell line.
- Additional viability data can be obtained from [DepMap's PRISM dataset](https://depmap.org/portal/download/).
- [Therapeutics Data Commons](https://github.com/mims-harvard/TDC) provides access to a number of compound databases as part of their [cheminformatics tasks](https://tdcommons.ai/benchmark/overview/). (In the same vein, [OpenProblems](https://openproblems.bio/) provides a [framework for tasks in single-cell](https://github.com/openproblems-bio/openproblems/tree/main/openproblems/tasks) which can also support perturbation modeling tasks in a more long term format than was previously seen in the [DREAM challenges](https://dreamchallenges.org/dream-7-nci-dream-drug-sensitivity-and-drug-synergy-challenge/).)

### Current approaches

We [maintain a list of perturbation-related tools at scrna-tools](https://www.scrna-tools.org/tools?sort=name&cats=Perturbations). Please consider further updating and tagging tools [there](https://github.com/scRNA-tools/scRNA-tools).

For the basis of the table in the article, see this [spreadsheet of a subset of perturbation models](https://docs.google.com/spreadsheets/d/1nqNg0DW1-Om7WtvRS20q-6b28usVRv5czOcxgj83Sgg/) includes more details.

### Datasets

Below, we curated a [table](https://raw.githubusercontent.com/theislab/sc-pert/main/data_table.csv) of perturbation datasets based on [Svensson *et al.* (2020)](https://doi.org/10.1093/database/baaa073).
