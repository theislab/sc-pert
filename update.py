import os
import glob
import pandas as pd

### data updates ###

# single-cell database
os.system('rm data.tsv')
os.system('rm personal.csv')
os.system('wget http://www.nxn.se/single-cell-studies/data.tsv')
os.system("wget --no-check-certificate -O personal.csv 'https://docs.google.com/spreadsheets/d/14awt-bCOnj4ca2uoKzuTNuKtUKXcoN82_-oGg2f1Ros/export?gid=1438063781&format=csv'")

# GDSC
#os.system('wget -P /gdsc/ ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC1_fitted_dose_response_25Feb20.xlsx')

path = 'datasets/'
processed_datasets = [file.split('.')[0] for file in glob.glob(f'{path}*.ipynb') if 'curation' not in file]
curated_datasets = [file.split('.')[0] for file in glob.glob(f'{path}*_curation.ipynb')]

personal_rec = pd.read_csv('personal.csv')
personal_rec['author_year'] = personal_rec['Author'] + '_' + personal_rec['Year'].astype(str)
dois = personal_rec.DOI.values

df = pd.read_csv('data.tsv', sep='\t')
df = df[df.DOI.isin(dois)]
df['Date'] = df['Date'].astype('object')  # ensure proper display
not_in_db = list(set(dois) - set(df.DOI.values))

# create placeholders for dois not in the sc studies DB
pdf = personal_rec[personal_rec.DOI.isin(not_in_db)]
df = df.append(pdf[pdf.columns.intersection(df.columns)])

# add additional columns of info
add_cols = personal_rec[['DOI', 'Treatment', '# perturbations', '# cell types', '# doses', '# timepoints', 'Author', 'Year']]
n_cols = 1 - add_cols.shape[1]  # default to negative
df = df.merge(
    add_cols,
    left_on='DOI',
    right_on='DOI'
)
df = df[list(df.columns[:5]) + \
    list(df.columns[n_cols:]) + list(df.columns[5:n_cols])]  # up to `Title`

# convert DOIs to links in markdown
links = []
for shorthand, link in df[['Shorthand', 'DOI']].values:
    s = f'[{shorthand}](https://doi.org/{link})'
    s = s.replace('et al', '*et al.*')
    links.append(s)
df['Shorthand'] = links

# add availability column with download links for curated datasets
links = []
for shorthand, author, year, link in df[['Shorthand', 'Author', 'Year', 'DOI']].values:
    s = ''
    base_nb_path = 'https://nbviewer.ipython.org/github/theislab/sc-pert/blob/main/datasets/'
    row = personal_rec[personal_rec.author_year == f'{author}_{year}']
    if not pd.isnull(row.Raw.values[0]):  # raw .h5ad
        s += f' [\\[raw h5ad\\]]({row.Raw.values[0]})'
    if not pd.isnull(row.Processed.values[0]):  # processed .h5ad
        s += f' [\\[processed h5ad\\]]({row.Processed.values[0]})'
    if f'{path}{author}_{year}_curation' in curated_datasets:  # curation notebook
        s += f' [\\[curation nb\\]]({base_nb_path}{author}_{year}_curation.ipynb)'
    if f'{path}{author}_{year}' in processed_datasets:  # processing notebook
        s += f' [\\[processing nb\\]]({base_nb_path}{author}_{year}.ipynb)'

    links.append(s)
df['.h5ad availability'] = links

# clean up
df = df.sort_values(by=['Treatment', 'Date'])
df = df.drop(['Authors', 'Journal', 'DOI', 'bioRxiv DOI', 'Author', 'Year', 'Date'], axis=1)

# rearrange columns
primary_cols = ['Shorthand', 'Title', '.h5ad availability', 'Treatment', '# perturbations', '# cell types', '# doses', '# timepoints']
df = df[primary_cols + \
    [c for c in df if c not in primary_cols]]

# write README
filenames = []
with open('README.md', 'w') as outfile:
    with open('readme_body.txt') as infile:
        outfile.write(infile.read())
        md = df.to_markdown(index=False, tablefmt='github', floatfmt='.8g')
        md = md.replace('| Title', '| Title'+'&nbsp;'*70)
        outfile.write(md)
        infile.close()
    outfile.close()

df.to_csv('data_table.csv')
