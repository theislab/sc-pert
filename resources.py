import os
import pandas as pd

def GDSC_response(cell_line, filename=None):
    """Downloads and returns the responsiveness ranking of compounds
    for a specific cell line from GDSC1 and GDSC2 as dataframes.

    Params
    ------
    cell_line : str
        Cell line identifier, e.g. A549.
    filename : str (default: None)
        The downloaded files will be saved as filename_GDSC1/2.csv.
    """
    if not os.path.exists('Cell_Lines_Details.xlsx'):
        os.system('wget ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Cell_Lines_Details.xlsx')
    annot = pd.read_excel('Cell_Lines_Details.xlsx')
    if cell_line not in annot['Sample Name'].values:
        raise ValueError(f'{cell_line} not present in GDSC')
    cell_id = annot[annot['Sample Name'] == cell_line]['COSMIC identifier'].values[0].astype('int').astype('str')
    if filename is None:
        filename = cell_line
    os.system(f"wget -O {filename}_GDSC1.csv 'https://www.cancerrxgene.org/api/cellline/download_zscore?id={cell_id}&screening_set=GDSC1&export=csv'")
    os.system(f"wget -O {filename}_GDSC2.csv 'https://www.cancerrxgene.org/api/cellline/download_zscore?id={cell_id}&screening_set=GDSC2&export=csv'")
    return pd.read_csv(f'{filename}_GDSC1.csv'), pd.read_csv(f'{filename}_GDSC2.csv')
