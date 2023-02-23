from io import StringIO

import streamlit as st
from serenipy.dtaselectfilter import results_from_df, to_dta_select_filter, DtaSelectFilterVersion, \
    from_dta_select_filter

import util

st.header('PaSER Diff! :bar_chart:')

st.write("""
This app is used to extract unique identifications from uploaded experiments. A resulting csv and DTASelect-filter file will be created
for each experiment, which will contain only peptides that were uniquely identified within that experiment, and no others. The Charge, modification, and protein
group options are used for set logic, and all charge/modifications will be included in the results files.""")

files = st.file_uploader(label='DTASelect-filter.txt files', accept_multiple_files=True, type='.txt')
use_charge, use_modifications, use_groups = util.peptide_config(disable_group=True)
group_by_peptide = st.checkbox(label='Group by peptide', value=True)

with st.expander('Custom Order'):
    orders, labels = util.get_file_order_and_labels(files)

if st.button('Run'):

    if len(set(orders)) != len(files):
        st.warning('Order must be unique!')
        st.stop()

    if len(set(labels)) != len(files):
        st.warning('Labels must be unique!')
        st.stop()

    orders, labels, files = zip(*sorted(zip(orders, labels, files)))

    versions, headers, results, tailers  = [], [], [], []
    for file in files:
        file_io = StringIO(file.getvalue().decode("utf-8"))
        version, header, dta_select_filter_results, tailer = from_dta_select_filter(file_io)
        versions.append(version)
        headers.append(header)
        results.append(dta_select_filter_results)
        tailers.append(tailer)

    df = util.get_peptides_and_protein_df(results)
    util.add_protein_groups(df, use_groups)
    util.add_peptide_groups(df, use_charge, use_modifications)

    with st.expander('Data'):
        st.dataframe(df)
        st.markdown(util.create_download_link(df.to_csv(index=False).encode('UTF-8'), f'combined.csv'),
                    unsafe_allow_html=True)

    protein_key_to_file_nums = {}
    peptide_key_to_file_nums = {}
    for protein_key, peptide_key, file_num in df[['protein_key', 'peptide_key', 'file_num']].values:
        protein_key_to_file_nums.setdefault(protein_key, set()).add(file_num)
        peptide_key_to_file_nums.setdefault(peptide_key, set()).add(file_num)

    KEY = 'peptide_key' if group_by_peptide is True else 'protein_key'

    st.header('Difference')
    for i in range(len(files)):
        st.subheader(labels[i])

        if group_by_peptide is True:
            df_diff = df[[len(peptide_key_to_file_nums[key]) == 1 and i in peptide_key_to_file_nums[key] for key in df['peptide_key']]]
        else:
            df_diff = df[[len(protein_key_to_file_nums[key]) == 1 and i in protein_key_to_file_nums[key] for key in df['protein_key']]]
        df_diff = df_diff[df_diff['file_num']==i]

        results = results_from_df(df_diff)
        results.sort(key=lambda x: x.protein_lines[0].sequence_coverage, reverse=True)
        filter_content = to_dta_select_filter(version=versions[i], h_lines=headers[i],
                                              dta_filter_results=results, end_lines=tailers[i])
        st.dataframe(df_diff)
        st.markdown(util.create_download_link(filter_content.encode('UTF-8'), f'{labels[i]}_diff.txt'),
                    unsafe_allow_html=True)
        st.markdown(util.create_download_link(df_diff.to_csv(index=False).encode('UTF-8'), f'{labels[i]}_diff.csv'),
                    unsafe_allow_html=True)

    st.header('Intersection')
    for i in range(len(files)):
        st.subheader(labels[i])
        if group_by_peptide is True:
            df_diff = df[[len(peptide_key_to_file_nums[key]) == len(files) and i in peptide_key_to_file_nums[key] for key in df['peptide_key']]]
        else:
            df_diff = df[[len(protein_key_to_file_nums[key]) == len(files) and i in protein_key_to_file_nums[key] for key in df['protein_key']]]

        df_diff = df_diff[df_diff['file_num']==i]
        results = results_from_df(df_diff)
        results.sort(key=lambda x: x.protein_lines[0].sequence_coverage, reverse=True)
        filter_content = to_dta_select_filter(version=versions[i], h_lines=headers[i],
                                              dta_filter_results=results, end_lines=tailers[i])
        st.dataframe(df_diff)
        st.markdown(util.create_download_link(filter_content.encode('UTF-8'), f'{labels[i]}_inter.txt'),
                    unsafe_allow_html=True)
        st.markdown(util.create_download_link(df_diff.to_csv(index=False).encode('UTF-8'), f'{labels[i]}_inter.csv'),
                    unsafe_allow_html=True)


