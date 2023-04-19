from io import StringIO

import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
from scipy.stats import stats
from serenipy.dtaselectfilter import from_dta_select_filter

import config
import util

st.header('PaSER Stats! :bar_chart:')

st.write("""
This app is used to generate venn diagrams for protein and peptide overlap between 2 or 3 experiments.
""")

with st.expander('Help'):
    st.markdown(config.PASER_VENN_HELP_MSG)

files = st.file_uploader(label='DTASelect-filter.txt files', accept_multiple_files=True, type='.txt')

with st.expander('Custom Order'):
    orders, labels = util.get_file_order_and_labels(files)

if st.button('Run'):

    if not len(files):
        st.warning('Upload files!')
        st.stop()

    orders, labels, files = zip(*sorted(zip(orders, labels, files)))
    data = []
    for i, file in enumerate(files):
        file_io = StringIO(file.getvalue().decode("utf-8"))
        _, _, results, _ = from_dta_select_filter(file_io)

        proteins = len([res.protein_lines[0].sequence_coverage for res in results])
        peptides = sum([len(res.peptide_lines) for res in results])
        data.append(
            {
                'name': labels[i],
                'order': str(i),
                'proteins': proteins,
                'peptides': peptides,
                'sequence_coverage': np.mean([res.protein_lines[0].sequence_coverage for res in results]),
                'sequence_coverage_norm': np.mean([res.protein_lines[0].sequence_coverage for res in results])*proteins,
                'spectrum_count': np.mean([res.protein_lines[0].spectrum_count for res in results]),
                'sequence_count': np.mean([res.protein_lines[0].sequence_count for res in results]),
                'nsaf': np.mean([res.protein_lines[0].nsaf for res in results]),
                'empai': np.mean([res.protein_lines[0].empai for res in results]),
                'x_corr': np.mean([line.x_corr for res in results for line in res.peptide_lines]),
                'delta_cn': np.mean([line.delta_cn for res in results for line in res.peptide_lines]),
                'conf': np.mean([line.conf for res in results for line in res.peptide_lines]),
                'sequence_coverage_sem': stats.sem([line.sequence_coverage for res in results for line in res.protein_lines]),
                'spectrum_count_sem': stats.sem([line.spectrum_count for res in results for line in res.protein_lines]),
                'sequence_count_sem': stats.sem([line.sequence_count for res in results for line in res.protein_lines]),
                'nsaf_sem': stats.sem([line.nsaf for res in results for line in res.protein_lines]),
                'empai_sem': stats.sem([line.empai for res in results for line in res.protein_lines]),
                'x_corr_sem': stats.sem([line.x_corr for res in results for line in res.peptide_lines]),
                'delta_cn_sem': stats.sem([line.delta_cn for res in results for line in res.peptide_lines]),
                'conf_sem': stats.sem([line.conf for res in results for line in res.peptide_lines])
            }
        )

    df = pd.DataFrame(data)

    fig = px.bar(df, x="order", y='sequence_coverage', text_auto=True, error_y='sequence_coverage_sem',
                 title='Average Protein Sequence Coverage',
                 labels={
                     "sequence_coverage": "Coverage %",
                     "order": ""
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y='sequence_coverage_norm', text_auto=True,
                 title='Average Normalized Protein Sequence Coverage',
                 labels={
                     "sequence_coverage_norm": "Normalized Coverage",
                     "order": ""
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y='spectrum_count', text_auto=True, error_y='spectrum_count_sem',
                 title='Average Protein Spectrum Count',
                 labels={
                     "spectrum_count": "Spectrum Count",
                     "order": ""
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y='sequence_count', text_auto=True, error_y='sequence_count_sem',
                 title='Average Protein Sequence Count',
                 labels={
                     "sequence_count": "Sequence Count",
                     "order": ""
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y='nsaf', text_auto=True, error_y='nsaf_sem',
                 title='Average Protein NSAF',
                 labels={
                     "nsaf": "NSAF",
                     "order": ""
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y='empai', text_auto=True, error_y='empai_sem',
                 title='Average Protein EMPAI',
                 labels={
                     "empai": "EMPAI",
                     "order": ""
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y='x_corr', text_auto=True, error_y='x_corr_sem',
                 title='Average Peptide XCORR',
                 labels={
                     "x_corr": "XCORR",
                     "order": ""
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y='delta_cn', text_auto=True, error_y='delta_cn_sem',
                 title='Average peptide Delta CN',
                 labels={
                     "delta_cn": "Delta CN",
                     "order": ""
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)

    fig = px.bar(df, x="order", y='conf', text_auto=True, error_y='conf_sem',
                 title='Average Peptide Confidence',
                 labels={
                     "conf'": "Confidence",
                     "order": ""
                 },
                 )
    fig.update_xaxes(tickangle=90,
                     tickmode='array',
                     tickvals=df['order'],
                     ticktext=df['name'])
    st.plotly_chart(fig)