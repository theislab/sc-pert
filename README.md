# sc-pert - Machine learning for perturbational single-cell omics

*This repository provides a community-maintained summary of models and datasets. It was initially curated for [(Cell Systems, 2021)](https://doi.org/10.1016/j.cels.2021.05.016).*

### External annotations

There are various resources for evaluation of single cell perturbation models. We discuss five tasks in the publication which can be supported by the following publicly available annotations:

- [GDSC](https://www.cancerrxgene.org/downloads/bulk_download) provides a collection of cell viability measurements for many compounds and cell lines. We provide a [code snippet](https://github.com/theislab/sc-pert/blob/main/resources.py#L4) to conveniently load GDSC-provided z-score compound response rankings per cell line.
- Additional viability data can be obtained from [DepMap's PRISM dataset](https://depmap.org/portal/download/).
- [Therapeutics Data Commons](https://github.com/mims-harvard/TDC) provides access to a number of compound databases as part of their [cheminformatics tasks](https://tdcommons.ai/benchmark/overview/). (In the same vein, [OpenProblems](https://openproblems.bio/) provides a [framework for tasks in single-cell](https://github.com/openproblems-bio/openproblems/tree/main/openproblems/tasks) which can also support perturbation modeling tasks in a more long term format than was previously seen in the [DREAM challenges](https://dreamchallenges.org/dream-7-nci-dream-drug-sensitivity-and-drug-synergy-challenge/).)
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) contains a comprehensive record of compounds ranging from experimental entities to non-proprietary small molecules. It is queryable via [PubChemPy](https://github.com/mcs07/PubChemPy).
- [DrugBank](https://go.drugbank.com/releases/latest) provides annotations for a relatively small number of small molecules in a standardized format.

### Current modeling approaches

We [maintain a list of perturbation-related tools at scrna-tools](https://www.scrna-tools.org/tools?sort=name&cats=Perturbations). Please consider further updating and tagging tools [there](https://github.com/scRNA-tools/scRNA-tools).

For the basis of the table in the article, see this [spreadsheet of a subset of perturbation models](https://docs.google.com/spreadsheets/d/1nqNg0DW1-Om7WtvRS20q-6b28usVRv5czOcxgj83Sgg/) which includes more details.

### Datasets

Below, we curated a [table](https://raw.githubusercontent.com/theislab/sc-pert/main/data_table.csv) of perturbation datasets based on [Svensson *et al.* (2020)](https://doi.org/10.1093/database/baaa073).

We also offer some datasets in a curated `.h5ad` format. These datasets have the following standardized fields in `.obs`:
* `perturbation_name` -- Human-readable ompound names (International non-proprietary naming where possible) for small molecules and gene names for genetic perturbations.
* `perturbation_type` -- `small molecule` or `genetic`
* `perturbation_value` -- A continuous covariate quantity, such as the dosage concentration or the number of hours since treatment.
* `perturbation_unit` -- Describes `perturbation_value`, such as `'ug'` or `'hrs'`.


| Shorthand                                                                  | Title&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                                                                                                                                                   | .h5ad availability                                                                                                                                                                                                                                     | Treatment       | # perturbations   | # cell types   | # doses   | # timepoints   |     Date | Reported cells total   | Organism     | Tissue               | Technique             | Data location   | Panel size   | Measurement   | Cell source                                           |   Disease | Contrasts             |   Developmental stage |   Number of reported cell types or clusters | Cell clustering   | Pseudotime   | RNA Velocity   | PCA   | tSNE   |   H5AD location | Isolation            | BC --> Cell ID _OR_ BC --> Cluster ID                              |   Number individuals |
|----------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------|-------------------|----------------|-----------|----------------|----------|------------------------|--------------|----------------------|-----------------------|-----------------|--------------|---------------|-------------------------------------------------------|-----------|-----------------------|-----------------------|---------------------------------------------|-------------------|--------------|----------------|-------|--------|-----------------|----------------------|--------------------------------------------------------------------|----------------------|
| [Jaitin *et al.* Science](https://doi.org/10.1126/science.1247651)         | Massively Parallel Single-Cell RNA-Seq for Marker-Free Decomposition of Tissues into Cell Types                                                         |                                                                                                                                                                                                                                                        | genetic targets | 8-22              | 1              | -         | 1              | 20140214 | 4,468                  | Mouse        | Spleen               | MARS-seq              | GSE54006        | nan          | RNA-seq       | CD11c+ enriched splenocytes                           |       nan | nan                   |                   nan |                                           9 | Yes               | No           | nan            | No    | No     |             nan | Sorting (FACS)       | nan                                                                |                  nan |
| [Adamson *et al.* Cell](https://doi.org/10.1016/j.cell.2016.11.048)        | A Multiplexed Single-Cell CRISPR Screening Platform Enables Systematic Dissection of the Unfolded Protein Response                                      |                                                                                                                                                                                                                                                        | genetic targets | 9-93 (sgRNA)      | 1              | -         | 1              | 20161215 | 86,000                 | Human        | Culture              | Perturb-seq           | GSE90546        | nan          | RNA-seq       | K562                                                  |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | Yes    |             nan | nan                  | nan                                                                |                  nan |
| [Dixit *et al.* Cell](https://doi.org/10.1016/j.cell.2016.11.038)          | Perturb-Seq: Dissecting Molecular Circuits with Scalable Single-Cell RNA Profiling of Pooled Genetic Screens                                            | [\[raw h5ad\]](https://ndownloader.figshare.com/files/34011689) [\[processed h5ad\]](https://ndownloader.figshare.com/files/34014608) [\[processing nb\]](https://nbviewer.ipython.org/github/theislab/sc-pert/blob/main/datasets/Dixit_2016.ipynb)    | genetic targets | 10,24             | 1              | -         | 1-2            | 20161215 | 200,000                | Human, Mouse | Culture              | Perturb-seq           | GSE90063        | nan          | RNA-seq       | BMDCs, K562                                           |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | No     |             nan | Nanodroplet dilution | nan                                                                |                  nan |
| [Datlinger *et al.* NMeth](https://doi.org/10.1038/nmeth.4177)             | Pooled CRISPR screening with single-cell transcriptome readout                                                                                          |                                                                                                                                                                                                                                                        | genetic targets | 3-29              | 1-2            | -         | 1              | 20170118 | 5,905                  | Human, Mouse | Culture              | CROP-seq              | GSE92872        | nan          | RNA-seq       | HEK293T, 3T3, Jurkat                                  |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | No     |             nan | nan                  | nan                                                                |                  nan |
| [Hill *et al.* NMethods](https://doi.org/10.1038/nmeth.4604)               | On the design of CRISPR-based single-cell molecular screens                                                                                             |                                                                                                                                                                                                                                                        | genetic targets | 32                | 1              | -         | 1              | 20180219 | 5,879                  | Human        | Culture              | CROP-seq              | GSE108699       | nan          | RNA-seq       | MCF10a cells                                          |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | https://github.com/shendurelab/single-cell-ko-screens#result-files |                  nan |
| [Gasperini *et al.* Cell](https://doi.org/10.1016/j.cell.2018.11.029)      | A Genome-wide Framework for Mapping Gene Regulation via Cellular Genetic Screens                                                                        |                                                                                                                                                                                                                                                        | genetic targets | 1119, 5779        | 1              | -         | 1              | 20190103 | 207,324                | Human        | Culture              | CROP-seq              | GSE120861       | nan          | RNA-seq       | K562 Cells                                            |       nan | CRISPR Screen         |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                  nan |
| [Jost *et al.* NBT](https://doi.org/10.1038/s41587-019-0387-5)             | Titrating gene expression using libraries of systematically attenuated CRISPR guide RNAs                                                                |                                                                                                                                                                                                                                                        | genetic targets | 25                | 2              | -         | 1              | 20200113 | 19,587                 | Human        | Culture              | Perturb-seq           | GSE132080       | nan          | RNA-seq       | K562 cells                                            |       nan | 25 gene screen        |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                  nan |
| [Schraivogel *et al.* NMethods](https://doi.org/10.1038/s41592-020-0837-5) | Targeted Perturb-seq enables genome-scale genetic screens in single cells                                                                               | [\[processing nb\]](https://nbviewer.ipython.org/github/theislab/sc-pert/blob/main/datasets/Schraivogel_2020.ipynb)                                                                                                                                    | genetic targets | 1778 (enhancers)  | 1              | -         | 1              | 20200601 | 231,667                | Human, Mouse | Bone marrow, Culture | TAP-seq               | GSE135497       | 1,000        | RNA-seq       | nan                                                   |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | Yes    |             nan | nan                  | nan                                                                |                  nan |
| [Ursu *et al.* bioRxiv](https://doi.org/10.1101/2020.11.16.383307)         | Massively parallel phenotyping of variant impact in cancer with Perturb-seq reveals a shift in the spectrum of cell states induced by somatic mutations |                                                                                                                                                                                                                                                        | genetic targets | 200               | 1              | -         | 1              | 20201118 | 162,314                | Human        | Lung                 | Perturb-seq           | nan             | nan          | RNA-seq       | nan                                                   |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                  nan |
| [Frangieh *et al.* NGenetics](https://doi.org/10.1038/s41588-021-00779-1)  | Multimodal pooled Perturb-CITE-seq screens in patient models define mechanisms of cancer immune evasion                                                 | [\[raw h5ad\]](https://ndownloader.figshare.com/files/34012565) [\[processed h5ad\]](https://ndownloader.figshare.com/files/34013717) [\[processing nb\]](https://nbviewer.ipython.org/github/theislab/sc-pert/blob/main/datasets/Frangieh_2021.ipynb) | genetic targets | 248               | 1              | -         | 1              | 20210301 | 218,331                | Human        | Culture              | Perturb-CITE-seq      | SCP1064         | nan          | RNA-seq       | nan                                                   |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                  nan |
| [Papalexi *et al.* NGenetics](https://doi.org/10.1038/s41588-021-00778-2)  | Characterizing the molecular regulation of inhibitory immune checkpoints with multimodal single-cell screens                                            |                                                                                                                                                                                                                                                        | genetic targets | nan               | nan            | nan       | nan            | 20210301 | 28,295                 | Human        | Culture              | CITE-seq & ECCITE-seq | GSE153056       | nan          | RNA-seq       | nan                                                   |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                  nan |
| [Datlinger *et al.* NMethods](https://doi.org/10.1038/s41592-021-01153-z)  | Ultra-high-throughput single-cell RNA sequencing and perturbation screening with combinatorial fluidic indexing                                         |                                                                                                                                                                                                                                                        | genetic targets | nan               | nan            | nan       | nan            | 20210531 | nan                    | Human, Mouse | nan                  | scifi-RNA-seq         | nan             | nan          | nan           | nan                                                   |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                  nan |
| [Norman *et al.* (2019)](https://doi.org/10.1126/science.aax4438)          | nan                                                                                                                                                     | [\[raw h5ad\]](https://ndownloader.figshare.com/files/34002548) [\[processed h5ad\]](https://ndownloader.figshare.com/files/34027562) [\[processing nb\]](https://nbviewer.ipython.org/github/theislab/sc-pert/blob/main/datasets/Norman_2019.ipynb)   | genetic targets | 278               | 1              | -         | 1              |      nan | nan                    | nan          | nan                  | CRISPRa               | nan             | nan          | RNA-seq       | induction of gene pair targets, single gene controls  |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                  nan |
| [Jin *et al.* (2020)](https://doi.org/10.1101/791525)                      | nan                                                                                                                                                     |                                                                                                                                                                                                                                                        | genetic targets | 35                | -              | -         | 1              |      nan | nan                    | nan          | nan                  | Perturb-seq           | nan             | nan          | RNA-seq       | in vivo mouse brain development                       |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                  nan |
| [Shin *et al.* SAdvances](https://doi.org/10.1126/sciadv.aav2249)          | Multiplexed single-cell RNA-seq via transient barcoding for simultaneous expression profiling of various drug perturbations                             |                                                                                                                                                                                                                                                        | small molecules | 45                | 2              | 1         | 1              | 20190516 | 3,091                  | Mouse, Human | Culture              | Drop-seq              | PRJNA493658     | nan          | RNA-seq       | HEK293T, NIIH3T3, A375, SW480, K562                   |       nan | 45 perturbations      |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                  nan |
| [Srivatsan *et al.* Science](https://doi.org/10.1126/science.aax6234)      | Massively multiplex chemical transcriptomics at single-cell resolution                                                                                  | [\[raw h5ad\]](https://ndownloader.figshare.com/files/33979517) [\[processing nb\]](https://nbviewer.ipython.org/github/theislab/sc-pert/blob/main/datasets/Srivatsan_2019.ipynb)                                                                      | small molecules | 188               | 3              | 4         | 2              | 20191206 | 650,000                | Human        | Culture              | sci-Plex              | GSE139944       | nan          | RNA-seq       | Cancer cell lines A549, K562, and MCF7                |       nan | 5,000 drug conditions |                   nan |                                           3 | Yes               | Yes          | No             | Yes   | No     |             nan | nan                  | nan                                                                |                  nan |
| [Zhao *et al.* bioRxiv](https://doi.org/10.1101/2020.04.22.056341)         | Deconvolution of Cell Type-Specific Drug Responses in Human Tumor Tissue with Single-Cell RNA-seq                                                       |                                                                                                                                                                                                                                                        | small molecules | 2,6               | 6,1            | -         | -              | 20200424 | 48,404                 | Human        | Brain, Tumor         | SCRB-seq (microwell)  | GSE148842       | nan          | RNA-seq       | nan                                                   |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                    6 |
| [McFarland *et al.* (2020)](https://doi.org/10.1101/868752)                | nan                                                                                                                                                     | [\[processing nb\]](https://nbviewer.ipython.org/github/theislab/sc-pert/blob/main/datasets/McFarland_2020.ipynb)                                                                                                                                      | small molecules | 1-13              | 24-99          | 1         | 1-5            |      nan | nan                    | nan          | nan                  | MIX-seq               | nan             | nan          | RNA-seq       | 4 small molecule experiments, one genetic             |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                  nan |
| [Chen *et al.* (2020)](https://doi.org/10.1038/s41592-019-0689-z)          | nan                                                                                                                                                     |                                                                                                                                                                                                                                                        | small molecules | 300               | 1              | 1         | 1              |      nan | nan                    | nan          | nan                  | CyTOF                 | nan             | nan          | protein       | breast  cancer  cells  undergoing  TGF-β-induced  EMT |       nan | nan                   |                   nan |                                         nan | nan               | nan          | nan            | nan   | nan    |             nan | nan                  | nan                                                                |                  nan |