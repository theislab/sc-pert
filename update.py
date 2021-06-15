import os
import pandas as pd

### data updates ###

# single-cell database
os.system('rm data.tsv')
os.system('wget http://www.nxn.se/single-cell-studies/data.tsv')

# GDSC
#os.system('wget -P /gdsc/ ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC1_fitted_dose_response_25Feb20.xlsx')

dois = [
### CRISPR ###
    '10.1016/j.cell.2016.11.038',  # Dixit 2016
    '10.1016/j.cell.2016.11.048',  # Adamson 2016
    '10.1126/science.1247651',     # Jaitin 2016
    '10.1038/nmeth.4177',          # Datlinger 2017
    '10.1038/nmeth.4604',          # Hill 2018
    '10.1016/j.cell.2018.11.029',  # Gasperini 2019
    # Norman 2019
    '10.1038/s41592-020-0837-5',   # Schraivogel 2020
    # Josh 2020
    '10.1101/2020.11.16.383307',   # Ursu 2020
    '10.1101/791525',              # Jin 2020
    '10.1038/s41588-021-00779-1',  # Frangieh 2021
### small molecules ### 
    '10.1126/science.aax6234',     # Srivatsan 2019
    '10.1126/sciadv.aav2249',      # Shin 2019
    '10.1101/868752',              # McFarland 2020
    # PhEMD dataset (Chen 2020) - cytobank
    '10.1101/2020.04.22.056341',   # Zhao 2020
    '10.1038/s41592-021-01153-z',  # Datlinger 2021
]

df = pd.read_csv('data.tsv', sep='\t')
df = df[df.DOI.isin(dois)]

# convert DOIs to links in markdown
links = []
for shorthand, link in df[['Shorthand', 'DOI']].values:
    s = f'[{shorthand}](doi.org/{link})'
    s = s.replace('et al', '*et al.*')
    links.append(s)
df['Shorthand'] = links
df = df.drop(['Authors', 'Journal', 'DOI', 'bioRxiv DOI'], axis=1)

# write README
filenames = []
with open('README.md', 'w') as outfile:
    with open('readme_body.txt') as infile:
        outfile.write(infile.read())
        md = df.to_markdown(index=False, tablefmt='github')
        md = md.replace('| Title', '| Title'+'&nbsp;'*100)
        outfile.write(md)
        infile.close()
    outfile.close()
