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
    st.markdown('''   
    
    **Peptide Uniqueness**
     
    * is defined by its **(sequence, charge)** pair. For example: (PEPTIDE, 2) != (PEPTIDE, 3) != (PEPT(15.2)IDE, 2) != (PEMMIDE, 2)
    
    **Data Frame Columns**
    
    * **name**: This element represents a label associated with the experiment.
    
    * **order**: This element represents the order of the experiment.
    
    * **unique_peptides**: This element represents the number of unique peptides in the experiment.
    
    * **duplicate_peptides**: This element represents the number of duplicate peptides in the experiment, i.e., the number of peptides that appear more than once in the experiment.
    
    * **total_peptides**: This element represents the total number of peptides in the experiment.
        
    * **new_unique_peptides**: This element represents the number of new unique peptides in the experiment that are not present in the previous experiments.
    
    * **new_peptides**: This element represents the number of new peptides in the experiment that are not present in the previous experiments.
    
    * **seen_unique_peptides**: This element represents the number of unique peptides in the experiment that are also present in the previous experiments.
    
    * **seen_peptides**: This element represents the number of peptides in the experiment that are also present in the previous experiments.
    
    * **unique_proteins**: This element represents the number of unique proteins in the experiment.
    
    * **duplicate_proteins**: This element represents the number of duplicate proteins in the data, i.e., the number of proteins that appear more than once in the experiment.
    
    * **total_proteins**: This element represents the total number of proteins in the experiment.
        
    * **new_proteins**: This element represents the number of new proteins in the experiment that are not present in the previous experiments.
    
    * **new_unique_proteins**: This element represents the number of new unique proteins in the experiment that are not present in the previous experiments.
    
    * **seen_proteins**: This element represents the number of proteins in the experiment that are also present in the previous experiments.
    
    * **seen_unique_proteins**: This element represents the number of unique proteins in the experiment that are also present in the previous experiments.
    ''')

files = st.file_uploader(label='DTASelect-filter.txt files', accept_multiple_files=True, type='.txt')

with st.expander('Custom Order'):

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

        data = {'name': label,
                'order': order,
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
