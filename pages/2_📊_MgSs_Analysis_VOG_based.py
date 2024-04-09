import streamlit as st
import pandas as pd
import numpy as np
import pickle
import plotly.graph_objects as go
from plotly.offline import plot

st.set_page_config(layout="centered")

######################################## ----- DATA IMPORTATION ----- ########################################

vog_clust = pd.read_csv("Data/mgss.clustering.parameters.vog.csv").rename(columns={"Unnamed: 0" : "species"}) # mgss clustering parameter
vog_mgss_coverage = pd.read_csv("Data/vog.mgSs.coverage.stats.csv").rename(columns={"Unnamed: 0" : "sub_species", "as.factor(sample_cluster)" : "sample_cluster"}) # mgss coverage stats
VOG_GeneProduct = pd.read_csv("Data/VOG_gene_product.csv")
# VOG_GeneProduct = pd.read_csv("/Users/amaros/Documents/VOG_gene_product.csv") #VIRGO2 GeneProduct with VOG annotation

# with open('/Users/amaros/Documents/reorder_dataframe.pkl', 'rb') as f:
with open('Data/reorder_dataframe.pkl', 'rb') as f:
    reorder_dataframe = pickle.load(f)

# with open("/Users/amaros/Documents/hover_dict.pkl", 'rb') as f:
with open('Data/hover.pkl', 'rb') as f:
        hover = pickle.load(f)

# with open('/Users/amaros/Documents/vog_mgss_pa.pkl', 'rb') as f:
with open('Data/vog_mgss_pa.pkl', 'rb') as f:
    vog_mgss_pa = pickle.load(f)

with open('Data/vog_clusters.pkl', 'rb') as f:
    vog_clusters = pickle.load(f)

# with open('/Users/amaros/Documents/gene_pa_count.pkl', 'rb') as f:
with open('Data/gene_pa_count.pkl', 'rb') as f:
    gene_count = pickle.load(f)


######################################## ----- DATA PROCESSING ----- ########################################

not_to_cluster = ["Alterileibacterium", "Anaerococcus", "Bacteroides", "Campylobacter",
                  "Corynebacterium", "Gardnerella", "Gulosibacter", "Lactobacillus", 
                  "Limosilactobacillus", "MultiGenera", "Porphyromonas", "Prevotella", "Streptococcus"]

vog_species = vog_mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]).unique()
species = [s for s in vog_species if s not in not_to_cluster]


######################################## ----- PAGE CONTENT ----- ########################################

st.title("Metagenomic Subspecies")

st.subheader("Visualizations")

option = st.selectbox("Species", species)

gene_count_df = gene_count[option]
Gene = pd.DataFrame(gene_count_df.set_index('Gene').loc[:,gene_count_df.columns[1]:].sum(axis = 1)).rename(columns = {0:"Number_of_samples"}).reset_index()
Gene['%_of_samples'] = Gene['Number_of_samples'].apply(lambda x : round((100 * x / (len(gene_count_df.columns) - 1)),2))

tab1, tab2, tab3 = st.tabs(["Species coverage", "Subspecies stats", "Presence Absence heatmap"])

with tab1:

        col1, col2 = st.columns(2)

        with col1 :
               st.image("Medias/vog_mgss_coverage_png/" + option + "_subspecies_coverage_boxplot.png")

        with col2 :
                st.image("Medias/vog_mgss_coverage_png/" + option + "_subspecies_coverage_by_NoVOG.png")

with tab2 :

        st.subheader("Clustering parameters")
        vog_df2 = vog_clust[vog_clust['species'] == option].drop('Number_of_unassigned', axis = 1)
        st.dataframe(vog_df2)

        st.subheader("Subspecies stats")
        st.dataframe(vog_mgss_coverage[vog_mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]) == option])

with tab3 :

        st.image("Medias/vog_heatmap_presence_absence/_" + option + "_heatmap_presence_absence.png")

        selection = st.radio("Select a filter", ["VOG cluster", "VOG ID"], index=None)

        if selection == "VOG cluster" :

                # Multiselection to choose one or multiple clusters of VOGs / ID
                vog_cluster_to_choose = st.multiselect("Select VOG cluster", [k for k in np.arange(1,11,1)])
                vog_cluster = [int(elt) for elt in vog_cluster_to_choose]

                vog_clusters_df = vog_clusters[option]
                vog_clusters_df['vog_cluster'] = vog_clusters_df['vog_cluster'].astype(int)
                vog_to_consider = vog_clusters_df[vog_clusters_df['vog_cluster'].isin(vog_cluster)]['vog']

                VOG = VOG_GeneProduct.loc[VOG_GeneProduct['VOG'].isin(vog_to_consider.values), :].reset_index(drop=True)
                st.dataframe(pd.merge(VOG, Gene, on='Gene', how='left'))

        if selection == "VOG ID":
                vog_id_to_choose = st.text_input("Select a VOG ID:")
                VOG = VOG_GeneProduct[VOG_GeneProduct['VOG'] == vog_id_to_choose]
                st.dataframe(pd.merge(VOG, Gene, on='Gene', how='left'))
        

        ### Interative Heatmap

        interactive_heatmap = st.button("Interactive heatmap")

        if interactive_heatmap :

                df_reorder = reorder_dataframe[option]
                
                # Create heatmap
                heatmap = go.Figure(data=go.Heatmap(
                z=df_reorder,
                x=df_reorder.columns,
                y=df_reorder.index,
                colorscale=[[0, 'antiquewhite'], [1, 'mediumblue']],
                showscale=False, 
                hovertext = hover[option],
                hovertemplate="VOG: %{y}<br>SampleID: %{x}<br>GeneProduct: %{hovertext}"

                ))

                # Update layout
                heatmap.update_layout(
                title= f'{option} Presence-Absence Heatmap',
                xaxis= dict(title='Samples'),
                yaxis=dict(title='VOG'),  
                )

                # Update x axis
                heatmap.update_xaxes(
                showticklabels=False
                )

                # Update y axis
                heatmap.update_yaxes(
                showticklabels=False
                )

                # Show plot
                # heatmap.show()
                plot(heatmap, filename='plot.html', auto_open=True)