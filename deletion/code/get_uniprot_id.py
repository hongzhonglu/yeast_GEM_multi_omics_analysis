import requests
from bs4 import BeautifulSoup
import time
from tqdm import tqdm, trange
import re
import pandas as pd


def get_name(uniprotid):
    headers = {
        'user-agent': "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/104.0.0.0 Safari/537.36",
    }
    url = 'https://rest.uniprot.org/uniprotkb/{}'.format(uniprotid)
    # url = 'https://rest.uniprot.org/uniprotkb/Q12499'

    response = requests.request("GET", url, headers=headers)
    rawdata = BeautifulSoup(response.text, 'html.parser')
    try:
        locus_names_r = re.findall(r'"value":"Y.*"', str(rawdata))
        locus_names = locus_names_r[0].split('"')[3]
    except IndexError:
        locus_names = 'nan'
    try:
        Standard_name_r = re.findall(r'"geneName":{"value":".*"}', str(rawdata))
        Standard_name = Standard_name_r[0].split('"')[5]
    except IndexError:
        Standard_name = 'nan'
    return locus_names, Standard_name

data = pd.read_csv('../data/yeast5k_noimpute_wide.csv')
data.index = data.loc[:, 'Protein.Group']
uni_id = list(data.index)
locus_names = []
standard_name = []
for i in trange(len(data.index)):
    locus, standard = get_name(uni_id[i])
    locus_names.append(locus)
    print(locus, standard)
    standard_name.append(standard)

data.index = locus_names

data.to_csv('../data/yeast5k_noimpute_wide.csv', index=True)