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
                'unique_peptides': len(peptide_set),
                'total_peptides': len(peptides),
                'duplicate_peptides': len(peptides) - len(peptide_set),
                'new_unique_peptides': len(peptide_set - global_peptides),
                'new_peptides': sum([p not in global_peptides for p in peptides]),
                'seen_unique_peptides': len(global_peptides.intersection(peptide_set)),
                'seen_peptides': sum([p in global_peptides for p in peptides]),
                'unique_proteins': len(protein_set),
                'duplicate_proteins': len(proteins) - len(protein_set),
                'total_proteins': len(proteins),
                'new_proteins': sum([p not in global_proteins for p in proteins]),
                'new_unique_proteins': len(protein_set - global_proteins),
                'seen_proteins': sum([p in global_proteins for p in proteins]),
                'seen_unique_proteins': len(global_proteins.intersection(protein_set))}

        data['unique_peptides_perc'] = round(data['unique_peptides'] / data['total_peptides'], 4) * 100
        data['duplicate_peptides_perc'] = round(data['duplicate_peptides'] / data['total_peptides'], 4) * 100
        data['new_peptides_perc'] = round(data['new_peptides'] / data['total_peptides'], 4) * 100
        data['seen_peptides_perc'] = round(data['seen_peptides'] / data['total_peptides'], 4) * 100
        data['new_proteins_perc'] = round(data['new_proteins'] / data['total_proteins'], 4) * 100
        data['seen_proteins_perc'] = round(data['seen_proteins'] / data['total_proteins'], 4) * 100
        global_proteins.update(protein_set)
        global_peptides.update(peptide_set)

        data['global_proteins'] = len(global_proteins)
        data['global_peptides'] = len(global_peptides)

        datas.append(data)

    df = pd.DataFrame(datas)
    st.dataframe(df)

    fig = px.bar(df, x="order", y=['unique_peptides', 'duplicate_peptides'], barmode="stack", text_auto=True,
                 title='Unique peptides vs duplicate peptides # (per experiment)')
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['unique_peptides_perc', 'duplicate_peptides_perc'], barmode="stack", text_auto=True,
                 title='Unique peptides vs duplicate peptides % (per experiment)')
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['new_peptides', 'seen_peptides'], barmode="stack", text_auto=True,
                 title='New peptides vs seen peptides # (per previous experiments)')
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['new_peptides_perc', 'seen_peptides_perc'], barmode="stack", text_auto=True,
                 title='New peptides vs seen peptides % (per previous experiments)')
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['new_proteins', 'seen_proteins'], barmode="stack", text_auto=True,
                 title='New proteins vs seen proteins # (per previous experiments)')
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['new_proteins_perc', 'seen_proteins_perc'], barmode="stack", text_auto=True,
                 title='New proteins vs seen proteins % (per previous experiments)')
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['global_proteins', 'global_peptides'], barmode="group", text_auto=True,
                 title='Global peptides and proteins (per previous experiments)')
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)
