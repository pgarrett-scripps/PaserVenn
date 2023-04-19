from io import StringIO

import streamlit as st

from matplotlib_venn import venn2
from matplotlib_venn import venn3
from matplotlib import pyplot as plt
from serenipy.dtaselectfilter import from_dta_select_filter

import config
import util

st.header('PaSER Venn! :bar_chart:')

st.write("""
This app is used to generate venn diagrams for protein and peptide overlap between 2 or 3 experiments.
""")

with st.expander('Help'):
    st.markdown(config.PASER_VENN_HELP_MSG)

files = st.file_uploader(label='DTASelect-filter.txt files', accept_multiple_files=True, type='.txt')
use_charge, use_modifications, use_groups = util.peptide_config()

with st.expander('Custom Order'):
    orders, labels = util.get_file_order_and_labels(files)

if st.button('Run'):

    if len(set(orders)) != len(files):
        st.warning('Order must be unique!')
        st.stop()

    if len(set(labels)) != len(files):
        st.warning('Labels must be unique!')
        st.stop()

    if len(files) != 2 and len(files) != 3:
        st.warning('Incorrect number of files: {len(files}. Please use only 2 or 3 files!')
        st.stop()

    orders, labels, files = zip(*sorted(zip(orders, labels, files)))

    results = []
    for file in files:
        file_io = StringIO(file.getvalue().decode("utf-8"))
        _, _, dta_select_filter_results, _ = from_dta_select_filter(file_io)
        results.append(dta_select_filter_results)

    df = util.get_peptides_and_protein_df(results)
    util.add_protein_groups(df, use_groups)
    util.add_peptide_groups(df, use_charge, use_modifications)

    with st.expander('Data'):
        #st.dataframe(df)
        st.markdown(util.create_download_link(df.to_csv(index=False).encode('UTF-8'), f'combined.csv'),
                    unsafe_allow_html=True)

    data = {}
    for i, grp in df.groupby('file_num'):
        data[i] = {'peptide': set(grp['peptide_key'].values),
                   'protein': set(grp['protein_key'].values)}

    protein_key_to_file_nums = {}
    peptide_key_to_file_nums = {}
    for protein_key, peptide_key, file_num in df[['protein_key', 'peptide_key', 'file_num']].values:
        protein_key_to_file_nums.setdefault(protein_key, set()).add(file_num)
        peptide_key_to_file_nums.setdefault(peptide_key, set()).add(file_num)

    df['protein_set'] = [protein_key_to_file_nums[protein_key] for protein_key in df['protein_key']]
    df['peptide_set'] = [peptide_key_to_file_nums[peptide_key] for peptide_key in df['peptide_key']]

    figure, axes = plt.subplots(2, 2)
    figure.tight_layout()

    protein_counts, peptide_counts = None, None
    total_peptides = None
    total_proteins = None
    if len(data) == 2:
        protein_sets = [data[0]['protein'],
                        data[1]['protein']]
        peptide_sets = [data[0]['peptide'],
                        data[1]['peptide']]
        protein_counts = [len(protein_set) for protein_set in protein_sets]
        peptide_counts = [len(peptide_set) for peptide_set in peptide_sets]

        v_protein = venn2(protein_sets, labels, ax=axes[0][0])
        v_peptide = venn2(peptide_sets, labels, ax=axes[1][0])

    elif len(data) == 3:
        protein_sets = [data[0]['protein'],
                        data[1]['protein'],
                        data[2]['protein']]
        peptide_sets = [data[0]['peptide'],
                        data[1]['peptide'],
                        data[2]['peptide']]

        protein_counts = [len(protein_set) for protein_set in protein_sets]
        peptide_counts = [len(peptide_set) for peptide_set in peptide_sets]

        v_protein = venn3(protein_sets, labels, ax=axes[0][0])
        v_peptide = venn3(peptide_sets, labels, ax=axes[1][0])

    else:
        st.warning('This should not have happened!')

    axes[0][1].barh(labels, protein_counts)
    axes[0][1].spines["top"].set_visible(False)
    axes[0][1].spines["right"].set_visible(False)
    axes[0][1].spines["bottom"].set_visible(False)
    axes[0][1].spines["left"].set_visible(False)

    axes[0][1].set_xticklabels([])
    axes[0][1].set_xticks([])

    for index, value in enumerate(protein_counts):
        axes[0][1].text(value, index, str(value))

    axes[1][1].barh(labels, peptide_counts)
    axes[1][1].spines["top"].set_visible(False)
    axes[1][1].spines["right"].set_visible(False)
    axes[1][1].spines["bottom"].set_visible(False)
    axes[1][1].spines["left"].set_visible(False)
    axes[1][1].set_xticklabels([])
    axes[1][1].set_xticks([])

    for index, value in enumerate(peptide_counts):
        axes[1][1].text(value, index, str(value))

    axes[0][1].title.set_text('Proteins')
    axes[1][1].title.set_text('Peptides')
    st.pyplot(fig=figure, clear_figure=None)

    st.markdown('---')
    st.subheader('Mapping')
    for i, file in enumerate(files):
        st.write(f'{labels[i]} -> {file.name}')

