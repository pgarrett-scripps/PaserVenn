from io import StringIO

import pandas as pd
import streamlit as st
from serenipy.dtaselectfilter import from_dta_select_filter
import plotly.express as px

st.header('PaSER Plot! :bar_chart:')

st.text("""
Generates plots to visualize peptide & protein IDs in sequential experiments.
""")

with st.expander('Help'):
    st.markdown('''Upload DTASelect-filter.txt files. Click Run.

    Unique peptides vs duplicate peptides (per experiment) - compares the number of unique (peptide, charge) pairs to duplicate (peptide, charge) pairs

    New peptides vs seen peptides (per previous experiments) - compares the number of new unique peptides found in sequential experiments. 
    
    New proteins vs seen proteins (per previous experiments) - ompares the number of new unique proteins found in sequential experiments. 
    ''')

files = st.file_uploader(label='DTASelect-filter.txt files', accept_multiple_files=True, type='.txt')

with st.expander('Custom Order'):
    labels = []
    orders = []
    for i, file in enumerate(files):
        name = file.name.split('.txt')[0]
        if len(file.name.split('_DTASelect-filter.txt')) > 1:
            name = file.name.split('_DTASelect-filter.txt')[0]
        elif len(file.name.split('DTASelect-filter.txt')) > 1:
            name = file.name.split('DTASelect-filter.txt')[0]

        st.caption(file.name)
        c1, c2 = st.columns(2)
        num = c1.number_input(label='Order', value=i + 1, key=f'num{file.name}')
        lab = c2.text_input(label='Label', value=name, key=f'lab{file.name}')

        orders.append(num)
        labels.append(lab)

if st.button('Run'):

    if len(set(orders)) != len(files):
        st.warning('Order must be unique!')
        st.stop()

    if len(set(labels)) != len(files):
        st.warning('Labels must be unique!')
        st.stop()

    orders, labels, files = zip(*sorted(zip(orders, labels, files)))

    peptide_lists, protein_lists = [], []
    for order, label, file in zip(orders, labels, files):
        peptides, proteins = [], []

        file_io = StringIO(file.getvalue().decode("utf-8"))
        _, _, dta_select_filter_results, _ = from_dta_select_filter(file_io)

        for res in dta_select_filter_results:
            for protein_line in res.protein_lines:
                proteins.append(protein_line.locus_name)

            for peptide_line in res.peptide_lines:
                sequence = peptide_line.sequence[2:-2]
                peptides.append((sequence, peptide_line.charge))

        peptide_lists.append(peptides)
        protein_lists.append(proteins)

    datas = []

    global_proteins = set()
    global_peptides = set()
    for peptides, proteins, label, order in zip(peptide_lists, protein_lists, labels, orders):
        peptide_set = set(peptides)
        protein_set = set(proteins)

        data = {}
        data['name'] = label
        data['order'] = order

        data['unique_peptides'] = len(peptide_set)
        data['total_peptides'] = len(peptides)
        data['duplicate_peptides'] = len(peptides) - len(peptide_set)
        data['new_peptides'] = len(peptide_set - global_peptides)
        data['seen_peptides'] = len(global_peptides.intersection(peptide_set))

        data['unique_proteins'] = len(protein_set)
        data['total_proteins'] = len(proteins)
        data['duplicate_proteins'] = len(proteins) - len(protein_set)
        data['new_proteins'] = len(protein_set - global_proteins)
        data['seen_proteins'] = len(global_proteins.intersection(protein_set))

        global_proteins.update(protein_set)
        global_peptides.update(peptide_set)

        data['total_unique_peptides'] = len(peptide_set)
        data['total_unique_proteins'] = len(protein_set)

        datas.append(data)

    df = pd.DataFrame(datas)
    st.dataframe(df)

    fig = px.bar(df, x="order", y=['unique_peptides', 'duplicate_peptides'], barmode="stack", text_auto=True,
                 title='Unique peptides vs duplicate peptides (per experiment)')
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['new_peptides', 'seen_peptides'], barmode="stack", text_auto=True,
                 title='New peptides vs seen peptides (per previous experiments)')
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y=['new_proteins', 'seen_proteins'], barmode="stack", text_auto=True,
                 title='New proteins vs seen proteins (per previous experiments)')
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)
