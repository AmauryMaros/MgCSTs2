import streamlit as st
import pandas as pd
from scipy.stats import kruskal
import numpy as np
from statsmodels.stats.multitest import multipletests
from itertools import combinations
import plotly.graph_objects as go
from sklearn.decomposition import PCA
import plotly.express as px

st.set_page_config(layout="wide")

######################################## ----- DATA IMPORTATION ----- ########################################

lsvf_raw = pd.read_csv("Data/LSVF_Metabolome_Raw.csv")
pressmat_raw = pd.read_csv("Data/PreSSMat_metabolomics_rawpeaks_100221.csv").rename(columns = {"UID":"sampleID"})
zapps_raw = pd.read_csv("Data/ZAPPS_metabolomics_rawpeaks_100221.csv").rename(columns = {"UID":"sampleID"})
mgCSTs_sort_vog = pd.read_csv("Data/vog_mgCST_sort_color.csv")
mgCSTs_sort_vog = mgCSTs_sort_vog[(mgCSTs_sort_vog['deepSplit'] == 4) & (mgCSTs_sort_vog['minClusterSize']==10)]
mgCSTs_samples_vog = pd.read_csv("Data/vog.mgCSTs.samples.df.csv")
mgCSTs_samples_vog = mgCSTs_samples_vog[(mgCSTs_samples_vog['deepSplit'] == 4) & (mgCSTs_samples_vog['minClusterSize']==10)]
zp_knn = pd.read_csv('Data/zapps_pressmat_clean.csv')
zp_min_imput = pd.read_csv('Data/zapps_pressmat_min_imputation.csv')
lsvf_knn = pd.read_csv("Data/lsvf_clean.csv")
lsvf_min_imput = pd.read_csv("Data/lsvf_min_imputation.csv")


######################################## ----- DATA PROCESSING ----- ########################################

lsvf_raw_merged = lsvf_raw.rename(columns = {"IGSbarcode":"sampleID"})
lsvf_raw_merged['sampleID'] = lsvf_raw_merged['sampleID'].astype(str).apply(lambda x : x + "_MG")
lsvf_raw_merged['project'] = "LSVF"
## PressMat Quality
rows_to_keep = ~pressmat_raw['sampleID'].str.contains("PBS")
pressmat_filtered = pressmat_raw[rows_to_keep].reset_index(drop = True)
pressmat_filtered['sampleID'] = pressmat_filtered['sampleID'].apply(lambda x : "MG_" + x)
pressmat_filtered.iloc[:,1:] = pressmat_filtered.iloc[:,1:].astype(float)
pressmat_filtered['project'] = "PressMat"
## ZAPPS Quality
rows_to_keep = ~zapps_raw['sampleID'].str.contains("PBS")
zapps_filtered = zapps_raw[rows_to_keep].reset_index(drop = True)
zapps_filtered['sampleID'] = zapps_filtered['sampleID'].apply(lambda x : "MG_" + x)
zapps_filtered.iloc[:,1:] = zapps_filtered.iloc[:,1:].astype(float)
zapps_filtered['project'] = "ZAPPS"
## Concatenate data
zapps_pressmat_raw = pd.concat([zapps_filtered, pressmat_filtered])
common_cols=[]
for cols in lsvf_raw_merged.columns:
     if cols in zapps_pressmat_raw.columns:
          common_cols.append(cols)

raw_merged = pd.merge(lsvf_raw_merged, zapps_pressmat_raw, on=common_cols, how='outer')
mgcsts_vog = mgCSTs_samples_vog[['sampleID', 'mgCST']]
raw_merged = pd.merge(raw_merged, mgcsts_vog, on='sampleID', how='inner')
raw_merged = raw_merged.dropna(subset='mgCST')

cols_to_remove = []
for i in raw_merged.columns:
    if raw_merged[i].isna().sum() == raw_merged.shape[0] :
        cols_to_remove.append(i)
raw_merged = raw_merged.drop(cols_to_remove, axis=1)

pca_data = raw_merged.drop(['sampleID','mgCST','project'],axis=1)
pca_data = pca_data.fillna(0.5*pca_data.min())
pca = PCA(n_components=3)
pca_data = pca.fit_transform(pca_data)
pca_data = pd.DataFrame(pca_data, columns=['PC1','PC2','PC3'])
pca_data['mgCST'] = raw_merged['mgCST']
pca_data['project'] = raw_merged['project']

df_var_norm = raw_merged.drop(['sampleID','mgCST','project'],axis=1)
df_var_norm = df_var_norm.fillna(0.5*df_var_norm.min())

# Calculate the variance of each feature
variances = df_var_norm.var()
# Normalize each feature by its variance
df_norm = df_var_norm.divide(variances.pow(0.5), axis=1)
pca_var = PCA(n_components=3)
pca_data_var = pca_var.fit_transform(df_norm)
pca_data_var = pd.DataFrame(pca_data_var, columns=['PC1','PC2','PC3'])
pca_data_var['mgCST'] = raw_merged['mgCST']
pca_data_var['project'] = raw_merged['project']

# Assign mgCST color for each sample
color_mgCST = mgCSTs_sort_vog[['mgCST', 'color_mgCST']].reset_index(drop = True)
color_mgCST['mgCST'] = color_mgCST['mgCST'].astype(str)
color_mgCST = dict(color_mgCST.values)

# Make MgCST columns as categorical column
pca_data['mgCST'] = pca_data['mgCST'].astype(str)
pca_data_var['mgCST'] = raw_merged['mgCST'].astype(str)

# Remove outlier on pca_data_var
pca_data_var = pca_data_var[pca_data_var['PC1'] <= 100]


######################################## ----- PAGE CONTENT ----- ########################################

st.subheader("Metabolome Analysis", divider = 'grey')

tab1, tab2 = st.tabs(['Non normalized data','Variance normalized data'])
with tab1 :
    col1, col2 = st.columns(2)
    with col1:

        fig = go.Figure()
        fig = px.scatter(pca_data, x='PC1', y='PC2', color='project',title="Project dependency - PCA representation")
        fig.update_layout(
        xaxis=dict(
            showgrid=True,  # Show grid lines
            gridcolor='rgba(0,0,0,0.1)',  # Grid line color
            gridwidth=1,  # Grid line width
        ))
        st.plotly_chart(fig)

        fig = go.Figure()
        fig = px.scatter(pca_data, x='PC2', y='PC3', color='project',title="Project dependency - PCA representation")
        fig.update_layout(
        xaxis=dict(
            showgrid=True,  # Show grid lines
            gridcolor='rgba(0,0,0,0.1)',  # Grid line color
            gridwidth=1,  # Grid line width
        ))
        st.plotly_chart(fig)

    with col2:

        fig = go.Figure()
        fig = px.scatter(pca_data, x='PC1', y='PC2', color='mgCST', color_discrete_map=color_mgCST, title="MgCSTs dependency - PCA representation")
        fig.update_layout(
        xaxis=dict(
            showgrid=True,  # Show grid lines
            gridcolor='rgba(0,0,0,0.1)',  # Grid line color
            gridwidth=1,  # Grid line width
        ))
        st.plotly_chart(fig)
        fig = go.Figure()
        fig = px.scatter(pca_data, x='PC2', y='PC3', color='mgCST', color_discrete_map=color_mgCST, title="MgCSTs dependency - PCA representation")
        fig.update_layout(
        xaxis=dict(
            showgrid=True,  # Show grid lines
            gridcolor='rgba(0,0,0,0.1)',  # Grid line color
            gridwidth=1,  # Grid line width
        ))
        st.plotly_chart(fig)

with tab2 :
    col1, col2 = st.columns(2)
    with col1:
        fig = px.scatter(pca_data_var, x='PC1', y='PC2', color='project',title="Project dependency - PCA representation")
        fig.update_layout(
        xaxis=dict(
            showgrid=True,  # Show grid lines
            gridcolor='rgba(0,0,0,0.1)',  # Grid line color
            gridwidth=1,  # Grid line width
        ))
        st.plotly_chart(fig)

        fig = px.scatter(pca_data_var, x='PC2', y='PC3', color='project',title="Project dependency - PCA representation")
        fig.update_layout(
        xaxis=dict(
            showgrid=True,  # Show grid lines
            gridcolor='rgba(0,0,0,0.1)',  # Grid line color
            gridwidth=1,  # Grid line width
        ))
        st.plotly_chart(fig)
    

    with col2:

        fig = go.Figure()
        fig = px.scatter(pca_data_var, x='PC1', y='PC2', color='mgCST', color_discrete_map=color_mgCST, title="MgCSTs dependency - PCA representation")
        fig.update_layout(
        xaxis=dict(
            showgrid=True,  # Show grid lines
            gridcolor='rgba(0,0,0,0.1)',  # Grid line color
            gridwidth=1))  # Grid line width
        st.plotly_chart(fig)
        fig = go.Figure()
        fig = px.scatter(pca_data_var, x='PC2', y='PC3', color='mgCST', color_discrete_map=color_mgCST, title="MgCSTs dependency - PCA representation")
        fig.update_layout(
        xaxis=dict(
            showgrid=True,  # Show grid lines
            gridcolor='rgba(0,0,0,0.1)',  # Grid line color
            gridwidth=1))  # Grid line width
        st.plotly_chart(fig)

st.subheader("MgCSTs distribution", divider='grey')

col1, col2 = st.columns(2)
with col1:
    value_counts = zp_knn['mgCST'].value_counts()
    fig = px.bar(x=value_counts.index, y=value_counts.values, title='MgCSTs distribution in Zapps PressMat')
    fig.update_layout(
    xaxis_title='MgCSTs',
    yaxis_title='Number of samples',
    xaxis=dict(
        tickmode='array',
        tickvals=value_counts.index,  # Set the tick values to the unique values in your DataFrame
        ticktext=value_counts.index,  # Set the tick labels to the unique values in your DataFrame
    )
)
    st.plotly_chart(fig)   
with col2:
    value_counts = lsvf_knn['mgCST'].value_counts()
    fig = px.bar(x=value_counts.index, y=value_counts.values, title='MgCSTs distribution in LSVF')
    fig.update_layout(
    xaxis_title='MgCSTs',
    yaxis_title='Number of samples',
    xaxis=dict(
        tickmode='array',
        tickvals=value_counts.index,  # Set the tick values to the unique values in your DataFrame
        ticktext=value_counts.index,  # Set the tick labels to the unique values in your DataFrame
    )
)
    st.plotly_chart(fig)

st.subheader("MgCSTs metabolome analysis", divider='grey')

col1, col2 =st.columns([3,2])

with col1:
    project_selection = st.selectbox(options=['LSVF','Zapps PressMat'], label="Select a project:")
    imputation_selection = st.selectbox(options=['kNN','Min value'], label="Select an imputation method:")
    mgcsts_selection = st.multiselect(options=range(1,27), label="Select at least 2 MgCSTs:")
    bullet_points = """
    * Double click on a label on the legend to isolate the corresponding plots
    * fold_change_A_vs_B : reference is A.
    * Log2(FC) is calculated by dividing A data by B data
    * Adjusted p-value has been estimated using Benjamini/Hochberg criterion
    """
    st.markdown(bullet_points)

with col2:
    st.write("MgCSTs dominant taxa:")
    mgCSTs_sort = mgCSTs_sort_vog.copy()
    mgCSTs_sort = mgCSTs_sort[(mgCSTs_sort['deepSplit'] == 4) & (mgCSTs_sort['minClusterSize'] == 10)].reset_index(drop=True).set_index('mgCST')
    st.dataframe(mgCSTs_sort['domTaxa'])

if (project_selection == "Zapps PressMat") and (imputation_selection == "kNN") :
     df = zp_knn.copy()
elif (project_selection == "Zapps PressMat") and (imputation_selection == "Min value") :
     df = zp_min_imput.copy()
elif (project_selection == "LSVF") and (imputation_selection == "Min value") :
     df = lsvf_min_imput.copy()
elif (project_selection == "LSVF") and (imputation_selection == "kNN") :
     df = lsvf_knn.copy()

if mgcsts_selection != []:
     
    df = df.drop('sampleID', axis=1)
    list_mgCSTs = sorted(mgcsts_selection)


    # List to store DataFrames for each mgCST value
    list_of_df = {}

    # Iterate through each unique mgCST value
    for mgCST_value in sorted(list_mgCSTs):
        # Extract DataFrame for the current mgCST value and drop 'mgCST' column
        df_mgCST = df[df['mgCST'] == mgCST_value].drop('mgCST', axis=1)
        # Append the DataFrame to the list
        list_of_df[mgCST_value] = df_mgCST


    # Function to calculate log2 fold change
    def log2_fold_change(group1, group2):
        return np.log2(group1.mean() / group2.mean())

    list_of_fc = ['fold_change_'+elt[0]+"_vs_"+elt[1] for elt in list(combinations([str(i) for i in list_mgCSTs],2))]

    log2_fc = {}
    for i in list_of_fc :
        a = int(i.split("_")[2])
        b = int(i.split("_")[4])
        df_a = df[df['mgCST'] == a].drop('mgCST', axis=1)
        df_b = df[df['mgCST'] == b].drop('mgCST', axis=1)

        log2_fc[i] = log2_fold_change(df_a, df_b)

    p_value_list = []
    adjusted_p_value_list = []
    # Iterate through each metabolite column
    cols = list_of_df[list(list_of_df.keys())[0]].columns
    for metabolite in cols:
        # Extract values of the metabolite column for each group
        values = [group[metabolite].values for group in list_of_df.values()]
        # Perform Kruskal-Wallis test
        h_stat, p_value = kruskal(*values)

        # # Perform multiple testing correction
        # _, adj_p_values_bh, _, _ = multipletests(p_value, method='fdr_bh')
        # _, adj_p_values_bonf, _, _ = multipletests(p_value, method='bonferroni')
        p_value_list.append(p_value)
        # adjusted_p_value_list.append(adj_p_values_bh[0])

    rejected, adjusted_p_values, _, _ = multipletests(p_value_list, method='fdr_bh')
    # adjusted_p_value_list.append(adjusted_p_values[0])

    result_stats = pd.DataFrame({"Metabolite":cols, "p_value":p_value_list, "adjusted_p_value" : adjusted_p_values})
    # result_stats = pd.DataFrame({"Metabolite":cols, "p_value":p_value_list})

    results_log2_fc = pd.DataFrame(log2_fc).reset_index().rename(columns={'index':'Metabolite'})

    results = pd.merge(result_stats, results_log2_fc, on='Metabolite', how='inner')


    fig = go.Figure()
            
    my_list = results.columns[3:]
    for i, key in enumerate(my_list):
            fig.add_trace(go.Scatter(x=results[key], 
                                    y=-np.log10(results['adjusted_p_value']),
                                    mode='markers',
                                    hoverinfo='text',
                                    text=results['Metabolite'],
                                    name=key))
            # Add vertical and horizontal lines
            fig.add_shape(type='line', x0=-1, x1=-1, y0=0, y1=max(-np.log10(results['adjusted_p_value'])), line=dict(color='green', width=2, dash='dash'))
            fig.add_shape(type='line', x0=1, x1=1, y0=0, y1=max(-np.log10(results['adjusted_p_value'])), line=dict(color='green', width=2, dash='dash'))
            fig.add_shape(type='line', x0=min(results[key]), x1=max(results[key]), y0=-np.log10(0.05), y1=-np.log10(0.05), line=dict(color='red', width=2, dash='dash'))


    fig.update_layout(
         xaxis=dict(
              showgrid=True,  # Show grid lines
              gridcolor='rgba(0,0,0,0.1)',  # Grid line color
              gridwidth=1,  # Grid line width
              ),
        legend=dict(
             orientation="h",  # Horizontal orientation
             yanchor="bottom",  # Anchor the legend to the bottom of the plot
             y=-0.4,  # Positioning of the legend below the plot
             xanchor="right",  # Anchor the legend to the right
             x=1  # Positioning of the legend to the right
             ),
             yaxis_title="-log10(adjusted_p_value)",
             xaxis_title="log2(FC)",
             width=1200,
             height=600)
    
    st.plotly_chart(fig, use_container_width=True)


    st.dataframe(results)