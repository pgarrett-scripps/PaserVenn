import base64
import re
from io import StringIO

import pandas as pd
import streamlit as st
from serenipy.dtaselectfilter import from_dta_select_filter, results_to_df


def create_download_link(val, filename):
    b64 = base64.b64encode(val)
    return f'<a href="data:application/octet-stream;base64,{b64.decode()}" download="{filename}">Download {filename}</a>'


def peptide_config(disable_charge=False, disable_mod=False, disable_group=False) -> (bool, bool, bool):
    use_charge = st.checkbox(label='Group by peptide charge', help='If False: (PEPTIDE +2 & PEPTIDE +3) == 1 unique peptides, '
                                                           'If True:  (PEPTIDE +2 & PEPTIDE +3) == 2 unique peptides',
                             value=True,
                             disabled=disable_charge)
    use_modifications = st.checkbox(label='Group by peptide modification',
                                    help='If False: (PEPTI(XXX)DE +2 & PEPTIDE +3) == 1 unique peptides, '
                                         'If True:  (PEPTI(XXX)DE +2 & PEPTIDE +3) == 2 unique peptides',
                                    value=True,
                                    disabled=disable_mod)
    use_groups = st.checkbox(label='Group by protein groups',
                             help='If False: all proteins in a group will be counted independently, '
                                  'If True: will count only the totla number of protein groups',
                             value=True,
                             disabled=disable_group)

    return use_charge, use_modifications, use_groups


def get_unmodified_peptide(peptide_sequence: str) -> str:
    pattern = re.compile(r'[^A-Z]')
    return pattern.sub('', peptide_sequence)


def get_file_order_and_labels(files):
    names = []
    for i, file in enumerate(files):
        name = file.name.split('.txt')[0]
        if len(file.name.split('_DTASelect-filter.txt')) > 1:
            name = file.name.split('_DTASelect-filter.txt')[0]
        elif len(file.name.split('DTASelect-filter.txt')) > 1:
            name = file.name.split('DTASelect-filter.txt')[0]

        names.append(name)

    if all([name.isalnum() for name in names]):
        pass

    labels = []
    orders = []
    for i, (file, name) in enumerate(zip(files, names)):
        st.caption(file.name)
        c1, c2 = st.columns(2)
        if name.isalnum():
            num = c1.number_input(label='Order', value=int(name), key=f'num{file.name}')
        else:
            num = c1.number_input(label='Order', value=i + 1, key=f'num{file.name}')
        lab = c2.text_input(label='Label', value=name, key=f'lab{file.name}')

        orders.append(num)
        labels.append(lab)

    return orders, labels


def get_peptides_and_protein_df(results):
    dfs = []
    for i, result in enumerate(results):
        tmp_df = results_to_df(result)
        tmp_df['file_num'] = i
        dfs.append(tmp_df)

    df = pd.concat(dfs)
    df['clean_sequence'] = df['sequence'].apply(lambda x: x[2:-2])
    df['unmod_sequence'] = df['clean_sequence'].apply(lambda x: get_unmodified_peptide(x))
    return df

def add_protein_groups(df, use_groups):
    seen = {}
    protein_keys = []
    for grp in df['protein_group' if use_groups is True else 'locus_name']:
        protein_keys.append(seen.setdefault(grp, len(protein_keys)))
    df['protein_key'] = protein_keys


def add_peptide_groups(df, use_charge, use_modifications):
    seen = {}
    peptide_keys = []
    for peptide, charge in df[['sequence' if use_modifications is True else 'unmod_sequence', 'charge']].values:
        if use_charge is True:
            peptide_keys.append(seen.setdefault((peptide, charge), len(seen)))
        else:
            peptide_keys.append(seen.setdefault(peptide, len(seen)))
    df['peptide_key'] = peptide_keys
