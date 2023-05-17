from io import StringIO

import pandas as pd
import streamlit as st
import plotly.express as px
from serenipy.dtaselectfilter import from_dta_select_filter

import config
import util

st.header('PaSER Plot! :bar_chart:')

st.text("""
This app is used to generate various plots to visualize peptide & protein IDs in sequential experiments.
""")

with st.expander('Help'):
    st.markdown(config.PASER_PLOT_HELP_MSG)

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
        st.dataframe(df)
        st.markdown(util.create_download_link(df.to_csv(index=False).encode('UTF-8'), f'combined.csv'),
                    unsafe_allow_html=True)

    datas = []
    global_proteins = set()
    global_peptides = set()
    for i in range(len(files)):
        tmp_df = df[df['file_num'] == i]
        peptides = list(tmp_df['peptide_key'].values)
        proteins = list(tmp_df['protein_key'].values)

        peptide_set = set(peptides)
        protein_set = set(proteins)

        data = {'name': labels[i],
                'order': str(i),
                'Unique Peptides': len(peptide_set),
                'Total Peptides': len(peptides),
                'Duplicate Peptides': len(peptides) - len(peptide_set),
                'New Unique Peptides': len(peptide_set - global_peptides),
                'New Peptides': sum([p not in global_peptides for p in peptides]),
                'Unique New Peptides': sum([p not in global_peptides for p in peptide_set]),
                'Seen Unique Peptides': len(global_peptides.intersection(peptide_set)),
                'Seen Peptides': sum([p in global_peptides for p in peptides]),
                'Unique Seen Peptides': sum([p in global_peptides for p in peptide_set]),
                'Unique Proteins': len(protein_set),
                'Duplicate Proteins': len(proteins) - len(protein_set),
                'Total Proteins': len(proteins),
                'New Proteins': sum([p not in global_proteins for p in protein_set]),
                'New Unique Proteins': len(protein_set - global_proteins),
                'Seen Proteins': sum([p in global_proteins for p in protein_set]),
                'Seen Unique Proteins': len(global_proteins.intersection(protein_set))}

        data['Unique Peptides Percent'] = round(data['Unique Peptides'] / data['Total Peptides'], 4) * 100
        data['Duplicate Peptides Percent'] = round(data['Duplicate Peptides'] / data['Total Peptides'], 4) * 100
        data['New Peptides Percent'] = round(data['New Peptides'] / data['Total Peptides'], 4) * 100
        data['Seen Peptides Percent'] = round(data['Seen Peptides'] / data['Total Peptides'], 4) * 100
        data['New Unique Peptides Percent'] = round(data['New Unique Peptides'] / data['Unique Peptides'], 4) * 100
        data['Seen Unique Peptides Percent'] = round(data['Seen Unique Peptides'] / data['Unique Peptides'], 4) * 100
        data['New Proteins Percent'] = round(data['New Proteins'] / data['Unique Proteins'], 4) * 100
        data['Seen Proteins Percent'] = round(data['Seen Proteins'] / data['Unique Proteins'], 4) * 100
        global_proteins.update(protein_set)
        global_peptides.update(peptide_set)

        data['Protein Counts'] = len(global_proteins)
        data['Peptide Counts'] = len(global_peptides)

        datas.append(data)

    df = pd.DataFrame(datas)
    st.dataframe(df)
    st.markdown(util.create_download_link(df.to_csv(index=False).encode('UTF-8'), f'paser_plot_results_{"_".join(labels)}.csv'),
                unsafe_allow_html=True)

    fig = px.bar(df, x="order", y=['Unique Peptides', 'Duplicate Peptides'], barmode="group", text_auto=True,
                 title='Number of unique peptides vs duplicate peptides in each experiment',
                 labels={
                     "unique_peptides": "Unique Peptide Count",
                     "duplicate_peptides": "Duplicate Peptide Count",
                     "order": "Experiment"
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['Unique Peptides Percent', 'Duplicate Peptides Percent'], barmode="stack", text_auto=True,
                 title='Percent of unique peptides vs duplicate peptides in each experiment',
                 labels={
                     "order": "Experiment"
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['New Peptides', 'Seen Peptides'], barmode="group", text_auto=True,
                 title='Number of new peptides vs seen previously seen peptides',
                 labels={
                     "order": "Experiment"
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['New Peptides Percent', 'Seen Peptides Percent'], barmode="stack", text_auto=True,
                 title='Percent of new peptides vs previously seen peptides',
                 labels={
                     "order": "Experiment"
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['New Unique Peptides', 'Seen Unique Peptides'], barmode="group", text_auto=True,
                 title='Number of new unique peptides vs seen previously seen unique peptides',
                 labels={
                     "order": "Experiment"
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['New Unique Peptides Percent', 'Seen Unique Peptides Percent'], barmode="stack", text_auto=True,
                 title='Percent of new unique peptides vs previously seen unique peptides',
                 labels={
                     "order": "Experiment"
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['New Proteins', 'Seen Proteins'], barmode="group", text_auto=True,
                 title='Number of new proteins vs previously seen proteins',
                 labels={
                     "order": "Experiment"
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['New Proteins Percent', 'Seen Proteins Percent'], barmode="stack", text_auto=True,
                 title='Percent of new proteins vs previously seen proteins',
                 labels={
                     "order": "Experiment"
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['Protein Counts', 'Peptide Counts'], barmode="group", text_auto=True,
                 title='Total number of peptides and proteins encountered in previous experiments',
                 labels={
                     "order": "Experiment"
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)
