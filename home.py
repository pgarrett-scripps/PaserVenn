import streamlit as st
from serenipy.dtaselectfilter import from_dta_select_filter

from matplotlib_venn import venn2
from matplotlib_venn import venn3
from matplotlib import pyplot as plt

st.header('PaSER Venn')

st.markdown('Upload 2-3 DTASelect-filter.txt files. Must have a unique name!')

files = st.file_uploader(label='DTASelect-filter.txt files', accept_multiple_files=True)

if st.button('Run'):

    if len(files) <= 1 or len(files) > 3:
        st.warning('Incorrect number of files: {len(files}. Please use only 2 or 3 files!')
        st.stop()

    data = {file.name: {'protein': set(), 'peptide':set()} for file in files}
    for file in files:
        version, head_lines, dta_select_filter_results, tail_lines = from_dta_select_filter(file.read().decode('utf-8'))

        for res in dta_select_filter_results:
            for protein_line in res.protein_lines:
                data[file.name]['protein'].add(protein_line.locus_name)
            for peptide_line in res.peptide_lines:
                data[file.name]['peptide'].add(peptide_line.sequence)

    print(data)
    figure, axes = plt.subplots(2, 1)
    if len(data) == 2:
        v_protein = venn2([data[list(data.keys())[0]]['protein'], data[list(data.keys())[1]]['protein']], ax=axes[0])
        v_peptide = venn2([data[list(data.keys())[0]]['peptide'], data[list(data.keys())[1]]['peptide']], ax=axes[1])

    if len(data) == 3:
        v_protein = venn3([data[list(data.keys())[0]]['protein'], data[list(data.keys())[1]]['protein'],
                          data[list(data.keys())[2]]['protein']], ax=axes[0])
        v_peptide = venn3([data[list(data.keys())[0]]['peptide'], data[list(data.keys())[1]]['peptide'],
                          data[list(data.keys())[2]]['peptide']], ax=axes[1])

    axes[0].title.set_text('Protein')
    axes[1].title.set_text('Peptide')
    st.pyplot(fig=figure, clear_figure=None)