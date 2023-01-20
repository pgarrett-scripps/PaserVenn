import re
from io import StringIO

import streamlit as st
from serenipy.dtaselectfilter import from_dta_select_filter

from matplotlib_venn import venn2
from matplotlib_venn import venn3
from matplotlib import pyplot as plt

st.header('PaserVenn')

st.text("""
This Streamlit app generates venn diagrams to visualize shared peptide/protein IDs from PaSER Experiments.
""")

with st.expander('Help'):
    st.markdown('''Upload 2-3 DTASelect-filter.txt files. Must have a unique name!
    
    Protein Counts - number of unique protein locus's
    
    Peptide Counts - number of unique peptide sequences
    ''')


files = st.file_uploader(label='DTASelect-filter.txt files', accept_multiple_files=True, type='.txt')
use_charge = st.checkbox(label='Consider charge', help='If False: (PEPTIDE +2 & PEPTIDE +3) == 1 unique peptides, '
                                                          'If True:  (PEPTIDE +2 & PEPTIDE +3) == 2 unique peptides')
use_modifications = st.checkbox(label='Consider modifications', help='If False: (PEPTI(XXX)DE +2 & PEPTIDE +3) == 1 unique peptides, '
                                                                  'If True:  (PEPTI(XXX)DE +2 & PEPTIDE +3) == 2 unique peptides')
use_groups = st.checkbox(label='Consider protein groups', help='If False: all proteins in a group will be counted independently, '
                                                                  'If True: will count only the totla number of protein groups')

def get_unmodified_peptide(peptide_sequence: str) -> str:
    pattern = re.compile(r'[^A-Z]')
    return pattern.sub('', peptide_sequence)


with st.expander('Custom Order'):
    labels = []
    order = []
    for i, file in enumerate(files):
        st.caption(file.name)
        c1, c2 = st.columns(2)
        num = c1.number_input(label='Order', value=i + 1, key=f'num{file.name}')
        lab = c2.text_input(label='Label', value=chr(65 + num - 1), key=f'lab{file.name}')

        order.append(num)
        labels.append(lab)

if st.button('Run'):

    if len(set(order)) != len(files):
        st.warning('Order must be unique!')
        st.stop()

    if len(set(labels)) != len(files):
        st.warning('Labels must be unique!')
        st.stop()

    order, labels, files = zip(*sorted(zip(order, labels, files)))

    if len(files) != 2 and len(files) != 3:
        st.warning('Incorrect number of files: {len(files}. Please use only 2 or 3 files!')
        st.stop()

    data = {file.name: {'protein': set(), 'peptide': set(), 'coverage': [], 'peptide_intensity': []} for file in files}
    for file in files:
        version, head_lines, dta_select_filter_results, tail_lines = from_dta_select_filter(StringIO(file.getvalue().decode("utf-8")))

        for res in dta_select_filter_results:
            for protein_line in res.protein_lines:
                data[file.name]['protein'].add(protein_line.locus_name)
                data[file.name]['coverage'].append(protein_line.sequence_coverage)

                if use_groups:
                    break

            for peptide_line in res.peptide_lines:

                sequence = peptide_line.sequence[2:-2]
                if not use_modifications:
                    sequence = get_unmodified_peptide(sequence)

                if use_charge:
                    data[file.name]['peptide'].add((sequence, peptide_line.charge))
                else:
                    data[file.name]['peptide'].add(sequence)

                data[file.name]['peptide_intensity'].append(peptide_line.total_intensity)

    shared_protein = set.intersection(*[data[file_name]['protein'] for file_name in data])
    shared_peptide = set.intersection(*[data[file_name]['peptide'] for file_name in data])
    unique_protein = set.union(*[data[file_name]['protein'] for file_name in data]) - set.intersection(
        *[data[file_name]['protein'] for file_name in data])
    unique_peptide = set.union(*[data[file_name]['peptide'] for file_name in data]) - set.intersection(
        *[data[file_name]['peptide'] for file_name in data])

    figure, axes = plt.subplots(2, 2)
    figure.tight_layout()

    protein_counts, peptide_counts = None, None
    if len(data) == 2:
        protein_sets = [data[list(data.keys())[0]]['protein'],
                        data[list(data.keys())[1]]['protein']]
        peptide_sets = [data[list(data.keys())[0]]['peptide'],
                        data[list(data.keys())[1]]['peptide']]
        protein_counts = [len(protein_set) for protein_set in protein_sets]
        peptide_counts = [len(peptide_set) for peptide_set in peptide_sets]

        v_protein = venn2(protein_sets, labels, ax=axes[0][0])
        v_peptide = venn2(peptide_sets, labels, ax=axes[1][0])

    elif len(data) == 3:
        protein_sets = [data[list(data.keys())[0]]['protein'],
                        data[list(data.keys())[1]]['protein'],
                        data[list(data.keys())[2]]['protein']]
        peptide_sets = [data[list(data.keys())[0]]['peptide'],
                        data[list(data.keys())[1]]['peptide'],
                        data[list(data.keys())[2]]['peptide']]
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

    axes[0][0].title.set_text('Shared Proteins')
    axes[1][0].title.set_text('Shared Peptides')
    axes[0][1].title.set_text('Protein Counts')
    axes[1][1].title.set_text('Peptide Counts')
    st.pyplot(fig=figure, clear_figure=None)

    st.markdown('---')
    st.subheader('Mapping')
    for i, file in enumerate(files):
        st.write(f'{labels[i]} -> {file.name}')
    st.markdown('---')
    st.subheader('Average Sequence Coverage')
    cols = st.columns(len(files))
    for i, file in enumerate(files):
        cols[i].metric(label=labels[i], value=str(round(sum(data[file.name]['coverage']) /
                                                        len(data[file.name]['coverage']), 2)) + "%")

    st.markdown('---')
    st.subheader('Average Peptide Intensity')
    cols = st.columns(len(files))
    for i, file in enumerate(files):
        cols[i].metric(label=labels[i], value=round(sum(data[file.name]['peptide_intensity']) /
                                                    len(data[file.name]['peptide_intensity']), 1))

    # figure = plt.figure()
    # plt.violinplot([data[file_name]['peptide_intensity'] for file_name in data])
    # plt.xticks(list(range(len(labels)+1)), [''] + labels)
    # plt.yscale('log')
    # st.pyplot(fig=figure, clear_figure=None)

    st.markdown('---')
    st.subheader('Stats')
    c1, c2, c3, c4 = st.columns(4)
    c1.metric(label='Shared proteins', value=len(shared_protein))
    c2.metric(label='Unique proteins', value=len(unique_protein))
    c3.metric(label='Shared peptides', value=len(shared_peptide))
    c4.metric(label='Unique peptides', value=len(unique_peptide))
