import streamlit as st
import pickle
import pandas as pd
from sklearn.decomposition import PCA
import plotly.express as px

st.set_page_config(layout="wide")


######################################## ----- CLUSTERING PARAMETERS SELECTION - SIDEBAR ----- ########################################

# Slider for parameters variation in the sidebar
st.sidebar.subheader("Parameters")
parameters = pd.read_csv("Data/mgCSTs_parameters_streamlit.csv")
value_ds = parameters['deepsplit'].values[0]
value_mcs = parameters['minClusterSize'].values[0]
deepsplit = st.sidebar.slider(label="deepSplit", min_value=0, max_value=4, value = value_ds)
minclustersize = st.sidebar.slider(label="minClusterSize", min_value=10, max_value=50, value = value_mcs)
parameters = pd.DataFrame({"minClusterSize" : [minclustersize], "deepsplit" : [deepsplit]})
parameters.to_csv("Data/mgCSTs_parameters_streamlit.csv", index = False)


######################################## ----- DATA IMPORTATION ----- ########################################

project = pd.read_csv("Data/VIRGO2_projects.csv")

@st.cache_data
def read_csv(url,minclustersize,deepsplit):
     df = pd.read_csv(url)
     df = df[(df['minClusterSize'] == minclustersize) & (df['deepSplit'] == deepsplit)]
     return df
mgcsts_samples_df = read_csv("Data/mgCSTs.samples.df.csv",minclustersize,deepsplit)
mgCSTs_sort = read_csv("Data/mgCST_sort_color.csv",minclustersize,deepsplit)

@st.cache_data
def read_pickle(path):
    with open(path, 'rb') as file:
        df = pickle.load(file)
    return df
# metabolomics = read_pickle("/Users/amaros/Documents/log_norm.pkl")
# pca_model =  read_pickle("Data/pca_model_no_phosphate_stearate.pkl")
metabolomics = read_pickle("Data/log_norm.pkl")
pca_model =  read_pickle("Data/pca_model.pkl")


@st.cache_data
def pca_model_data(minclustersize, deepsplit):
    return pca_model[minclustersize][deepsplit][0]
principal_components = pca_model_data(minclustersize, deepsplit)['principal_components']
explained_variance = pca_model_data(minclustersize, deepsplit)['explained_var_ratio']
sampleID  = pca_model_data(minclustersize, deepsplit)['sampleID']
mgCST = pca_model_data(minclustersize, deepsplit)['mgCST']


######################################## ----- DATA PROCESSING ----- ########################################

# Remove phosphate and stearate data from metabolomic columns
cols_to_drop = []
for i in metabolomics.columns[1:]:
    if ":phosphate" in i :
        cols_to_drop.append(i)
    elif ":stearate" in i :
        cols_to_drop.append(i)
metabolomics = metabolomics.drop(cols_to_drop, axis=1)

# Work on a copy of imported dataframe
mgcsts_samples = mgcsts_samples_df
mgcsts = mgCSTs_sort

# Create dataframe to show differences in number of samples between all/metabolome data
mgcsts = mgcsts.reset_index(drop = True)
count_sample = []
for element in mgcsts['dtc'].values :
    count_sample.append(mgcsts_samples.groupby(['dtc']).count()['sampleID'][element])
mgcsts['count_sample'] = count_sample

# Assign mgCST color for each sample
color_mgCST = mgcsts[['mgCST', 'color_mgCST']].reset_index(drop = True)
color_mgCST = color_mgCST[color_mgCST['mgCST'].isin(mgCST)]

# Create a dataframe with PCA datas for visualization
principal_components = pd.DataFrame(data=principal_components, columns = ['PC1','PC2','PC3'])
pca_df = pd.concat([sampleID, principal_components,mgCST], axis=1)
pca_df = pd.merge(pca_df, project, on='sampleID', how='inner') 
pca_df = pd.merge(pca_df, color_mgCST, on='mgCST', how='inner')      
pca_df = pca_df.sort_values(by='mgCST', ascending=True)
pca_df = pca_df.sort_values('mgCST', ascending=True).reset_index(drop=True)
pca_df['mgCST'] = pca_df['mgCST'].astype(str)   # to be interpreted as categorical column for plotly



######################################## ----- PAGE CONTENT ----- ########################################

st.header("Principal Component Analysis", divider='gray')

st.container()

col1, col2 = st.columns(2)

# PCA plot colored by MgCSTs
with col1 :

    fig = px.scatter(
        pca_df,
        x='PC1',
        y='PC2',
        color='mgCST',
        color_discrete_sequence=color_mgCST['color_mgCST'].values,
        title="PCA - colored by mgCSTs",
        labels = {'PC1' : "PC1 : "+str(round(explained_variance[0]*100,2))+"%",
                  'PC2' : "PC2 : "+str(round(explained_variance[1]*100,2)) + "%"})

    fig.add_annotation(
        dict(text=f"MinClusterSize = {minclustersize} <br> DeepSplit = {deepsplit}",
             xref='paper',
             yref='paper',
             x=0, y=1.05,
             showarrow=False,
             font=dict(size=12)
            ))
    
    st.plotly_chart(fig, use_container_width=True)

# PCA plot colored by project
with  col2 :

    fig = px.scatter(
        pca_df, 
        x='PC1', 
        y='PC2', 
        color='Project',
        title="PCA - colored by projects",
        labels = {'PC1' : "PC1 : "+str(round(explained_variance[0]*100,2))+"%",
                  'PC2' : "PC2 : "+str(round(explained_variance[1]*100,2)) + "%"})
    
    fig.add_annotation(
        dict(
            text=f"MinClusterSize = {minclustersize} <br> DeepSplit = {deepsplit}",
            xref='paper', yref='paper',
            x=0, y=1.05,
            showarrow=False,
            font=dict(size=12)
            ))
    
    st.plotly_chart(fig, use_container_width=True)


st.container()

col1, col2 = st.columns(2)

# Superposition of 2 tables to compare number of samples in metabolome data and clustered data (we don't have metabolome data for all samples)
with col1 :

    st.subheader("Number of samples in metabolomic data VS original data", divider='gray')
    st.write("Metabolomic")
    st.dataframe(pd.DataFrame(mgCST.value_counts().sort_index().rename(index = "count_sample")).transpose())
    st.write("Original")
    st.dataframe(mgcsts[['mgCST','count_sample']].set_index('mgCST').transpose())


# 3D PCA plot
st.container()
with col2 :

    st.subheader("3D Representation", divider='gray')

    tab1, tab2 = st.tabs(["MgCST", "Project"])

    # Plot colored by MgCST
    with tab1 :
        fig = px.scatter_3d(
            pca_df, 
            x='PC1', 
            y='PC2', 
            z='PC3', 
            color='mgCST',
            color_discrete_sequence=color_mgCST['color_mgCST'].values,
            # title="PCA - colored by mgCSTs",
            labels = {'PC1' : "PC1 : "+str(round(explained_variance[0]*100,2))+"%",
                      'PC2' : "PC2 : "+str(round(explained_variance[1]*100,2)) + "%",
                      'PC3' : "PC3 : "+str(round(explained_variance[2]*100,2))+"%"})
        
        fig.update_layout(legend= {'itemsizing': 'constant'})

        fig.update_traces(marker=dict(size=2, line=dict(width=1)),
                          selector=dict(mode='markers'))
        
        st.plotly_chart(fig, use_container_width=True)

    # Plot colored by Project
    with tab2 :
        fig = px.scatter_3d(
            pca_df, 
            x='PC1', 
            y='PC2', 
            z='PC3', 
            color='Project',
            # title="PCA - colored by projects",
            labels = {'PC1' : "PC1 : "+str(round(explained_variance[0]*100,2))+"%",
                      'PC2' : "PC2 : "+str(round(explained_variance[1]*100,2)) + "%",
                      'PC3' : "PC3 : "+str(round(explained_variance[2]*100,2))+"%"})
        
        fig.update_layout(legend= {'itemsizing': 'constant'})

        fig.update_traces(marker=dict(size=2, line=dict(width=1)),
                          selector=dict(mode='markers'))
        
        st.plotly_chart(fig, use_container_width=True)



st.header("Comparison between 2 groups", divider='gray')

# Create a double-slider to select 2 groups of data (it will be used to create a PCA)
st.container()
borne = mgcsts['mgCST'].nunique()
col1, col2 = st.columns(2)
with col1 :
    mgCST1 = st.slider("GroupA", 1, borne,(1,borne), key='mgCST1')    # this return a tuple (value1, value2)
with col2 :
    mgCST2 = st.slider("GroupB", 1, borne,(1,borne), key='mgCST2')

# Function to perform PCA on subset of data from double-slider selection
@st.cache_data
def data_pca(m1,m2):
    df1 = mgcsts_samples[mgcsts_samples['mgCST'].isin(range(m1[0],m1[1]+1))]
    df1.loc[:,'label'] = "GroupA"
    df2 = mgcsts_samples[mgcsts_samples['mgCST'].isin(range(m2[0],m2[1]+1))]
    df2.loc[:,'label'] = "GroupB"
    df = pd.concat([df1,df2], axis = 0)
    data_1 = pd.merge(df, metabolomics, on = "sampleID", how = "inner")
    return data_1

data1 = data_pca(mgCST1, mgCST2)


mgCSTs = data1['mgCST']
groups = data1['label']
data1 = data1.drop(['dtc','domTaxa','relabund','minClusterSize','deepSplit', 'sampleID', 'mgCST', 'label'], axis = 1)
new_colors = color_mgCST[color_mgCST['mgCST'].isin(mgCSTs)]


run_pca = st.button('Run PCA', key='run_pca')

if run_pca :

    # Create and fit the PCA model
    pca = PCA(n_components=6)
    principal_components = pca.fit_transform(data1)
    explained_variance = pca.explained_variance_ratio_
    components_compo = pca.components_

    # Create a dataframe with new PCA data
    pca_df = pd.DataFrame(data = principal_components, columns=['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6'])
    pca_df = pd.concat([pca_df, mgCSTs], axis=1)
    pca_df = pca_df.sort_values(by='mgCST', ascending=True)
    new_legend = new_colors['color_mgCST']#.apply(lambda x : mcolors.to_rgba(x)).values
    pca_df['mgCST'] = pca_df['mgCST'].astype(str)
    
    # Display PCA representation
    @st.cache_data
    def display_pca(df, pc0, pc1, a ,b):

        fig = px.scatter(
            df, 
            x=pc0, 
            y=pc1, 
            color='mgCST', 
            symbol='mgCST',
            color_discrete_sequence=new_legend.values,
            title=pc0+" - "+pc1+" representation",
            labels = {pc0 : pc0 + " : " + str(round(explained_variance[a]*100,2))+"%",
                      pc1 : pc1 + " : " + str(round(explained_variance[b]*100,2)) + "%"})
        
        fig.add_annotation(
            dict(
                text=f"MinClusterSize = {minclustersize} <br> DeepSplit = {deepsplit}",
                xref='paper', yref='paper',
                x=0, y=1.05,
                showarrow=False,
                font=dict(size=12)
                ))
        return fig
       
    cols = metabolomics.columns
    

    # Display 3 PCA representation side by side
    st.subheader("PCA Visualization")

    col1, col2, col3 = st.columns(3)

    with col1 :
        st.plotly_chart(display_pca(pca_df, 'PC1', 'PC2', 0, 1), use_container_width=True)
    with col2 :
        st.plotly_chart(display_pca(pca_df, 'PC3', 'PC4', 2, 3), use_container_width=True)
    with col3 :
        st.plotly_chart(display_pca(pca_df, 'PC5', 'PC6', 4, 5), use_container_width=True)

    # Features importances representation
    st.subheader("Features importances (explained variance)")

    n = 5     # Number of features to consider

    # Principal Components Loadings (original features composition)
    features_names = data1.columns
    principal_components_loadings = pd.DataFrame(pca.components_, columns=features_names)
    
    def get_features(pc, n) :

        '''For a given principal component, this function returns the n most contributing features (positive & negative correlation) from principal_components_loadings dataset.
        Results are stored in a dataframe composed of feature name, explained_variance associated, principal component desired and feature rank'''

        df_pos = pd.DataFrame(principal_components_loadings.iloc[pc,:].sort_values(ascending=False)[:n]).reset_index().rename(columns={'index':'Features', pc:'Explained_variance'})
        df_pos['PCs'] = 'PC'+str(pc+1)
        df_pos['Feature_rank'] = [i+1 for i in range(n)]

        df_neg = pd.DataFrame(principal_components_loadings.iloc[pc,:].sort_values(ascending=False)[-n:]).reset_index().rename(columns={'index':'Features', pc:'Explained_variance'})
        df_neg['PCs'] = 'PC'+str(pc+1)
        df_neg['Feature_rank'] = sorted([i+1 for i in range(n)],reverse=True)

        return df_pos, df_neg
    
    # Store positive & negative correlation DataFrame into a list
    top_pos = [get_features(i, 5)[0] for i in range(6)]
    top_neg = [get_features(i, 5)[1] for i in range(6)]

    # Concatenate all element of each list one after the other (vertical axis)
    plotly_df_pos = pd.DataFrame()
    for i in top_pos :
        plotly_df_pos = pd.concat([plotly_df_pos,i], axis=0)

    plotly_df_neg = pd.DataFrame()
    for i in top_neg :
        plotly_df_neg = pd.concat([plotly_df_neg, i], axis = 0)

    # Display the barplot for features loadings
    col1, col2 = st.columns(2)

    with col1 :

        plotly_df_pos['Feature_rank'] = plotly_df_pos['Feature_rank'].astype(str)

        fig = px.bar(
            plotly_df_pos, 
            x='PCs', 
            y='Explained_variance', 
            color='Feature_rank', 
            text='Features', 
            barmode='stack',
            labels={'PCs':'Principal Component'}, 
            title='5 most contributing features - Positive correlation')
        
        fig.add_annotation(
            dict(
                text=f"MinClusterSize = {minclustersize} <br> DeepSplit = {deepsplit}",
                xref='paper', yref='paper',
                x=0, y=1.05,
                showarrow=False,
                font=dict(size=12)
                ))
        
        st.plotly_chart(fig, theme='streamlit', use_container_width=True)

    with col2 :

        plotly_df_neg['Feature_rank'] = plotly_df_neg['Feature_rank'].astype(str)

        fig = px.bar(
            plotly_df_neg, 
            x='PCs', 
            y = 'Explained_variance', 
            color='Feature_rank',
            text='Features',
            barmode='stack',
            labels={'PCs':'Principal Component'}, 
            title='5 most contributing features - Negative correlation')

        fig.add_annotation(
            dict(
                text=f"MinClusterSize = {minclustersize} <br> DeepSplit = {deepsplit}",
                xref='paper', yref='paper',
                x=0, y=1.05,
                showarrow=False,
                font=dict(size=12)
                ))
        
        st.plotly_chart(fig, theme='streamlit', use_container_width=True)
