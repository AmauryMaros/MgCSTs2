import streamlit as st
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import pickle

st.set_page_config(layout="centered")

# VOG based data
vog_clust = pd.read_csv("Data/mgss.clustering.parameters.vog.csv").rename(columns={"Unnamed: 0" : "species"})
vog_mgss_coverage = pd.read_csv("Data/vog.mgSs.coverage.stats.csv").rename(columns={"Unnamed: 0" : "sub_species", "as.factor(sample_cluster)" : "sample_cluster"})
vog_species = vog_mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]).unique()

VOG_GeneProduct = pd.read_csv("Data/VOG_gene_product.csv")

not_to_cluster = ["Alterileibacterium", "Anaerococcus", "Bacteroides", "Campylobacter", "Corynebacterium", "Gardnerella", "Gulosibacter", "Lactobacillus", "Limosilactobacillus", "MultiGenera", "Porphyromonas", "Prevotella", "Streptococcus"]
species = [s for s in vog_species if s not in not_to_cluster]

with open('Data/vog_clusters.pkl', 'rb') as f:
    vog_clusters = pickle.load(f)

with open('Data/vog_mgss_pa.pkl', 'rb') as f:
    vog_mgss_pa = pickle.load(f)

with open('Data/samples_clusters_vog.pkl', 'rb') as f:
    samples_clusters_vog = pickle.load(f)

st.title("Metagenomic Subspecies")

st.subheader("Visualizations")

# a = len(species)
# st.write(a)
option = st.selectbox("Species", species)
tab1, tab2, tab3 = st.tabs(["Species coverage", "Subspecies stats", "Presence Absence heatmap"])

with tab1:
        col1, col2 = st.columns(2)
        with col1 :
        #     fig1 = Image.open("Medias/mgss_coverage_png/" + option + "_subspecies_coverage_boxplot.png") 
        #     st.image(fig1)
                st.image("Medias/vog_mgss_coverage_png/" + option + "_subspecies_coverage_boxplot.png")

        with col2 :
        #     fig2 = Image.open("Medias/mgss_coverage_png/" + option + "_subspecies_coverage_by_NoGenes.png")
        #     st.image(fig2)
                st.image("Medias/vog_mgss_coverage_png/" + option + "_subspecies_coverage_by_NoVOG.png")

with tab2 :
        st.subheader("Clustering parameters")
        vog_df2 = vog_clust[vog_clust['species'] == option].drop('Number_of_unassigned', axis = 1)
        st.dataframe(vog_df2)

        st.subheader("Subspecies stats")
        st.dataframe(vog_mgss_coverage[vog_mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]) == option])

with tab3 :
        st.image("Medias/vog_heatmap_presence_absence/_" + option + "_heatmap_presence_absence.png")

        # Multiselection to choose one or multiple clusters of VOGs
        vog_cluster_to_choose = st.multiselect("Select VOG cluster", [k for k in np.arange(1,11,1)])
        vog_cluster = [int(elt) for elt in vog_cluster_to_choose]
        
        vog_clusters_df = vog_clusters[option]
        vog_clusters_df['vog_cluster'] = vog_clusters_df['vog_cluster'].astype(int)
        vog_to_consider = vog_clusters_df[vog_clusters_df['vog_cluster'].isin(vog_cluster)]['vog']

        # # Multiselection to choose one or multiple clusters of samples
        # samples_cluster_to_choose = st.multiselect("Select samples cluster", [k for k in np.arange(1,11,1)])
        # samples_cluster = [int(elt) for elt in samples_cluster_to_choose]

        # samples_ids = samples_clusters_vog[option]
        # samples_ids['sample_cluster'] = samples_ids['sample_cluster'].astype(int)
        # samples_ids = samples_ids[samples_ids['sample_cluster'].isin(samples_cluster)]
        # samples_to_consider = samples_ids['sampleID']

        # # Presence absence selection
        # presence_absence = st.multiselect("Choose presence or absence:", [0,1])
        # # Filter the VOG presence absence table
        # vog_mgss_pa = vog_mgss_pa[option]
        # columns_to_keep = ['VOG'] + samples_to_consider.to_list()
        # vog_mgss_pa_filter = vog_mgss_pa[vog_mgss_pa['VOG'].isin(vog_to_consider)][columns_to_keep]
        # vog_mgss_pa_filter
        # # vog_to_consider = vog_mgss_pa_filter['VOG']

        # Get the dataframe from dictionary pickle for the considered mgss
        vog_clusters_df = vog_clusters[option]
        vog_clusters_df['vog_cluster'] = vog_clusters_df['vog_cluster'].astype(int)
        # Get the VOG Ids based on VOG_cluster selection in Streamlit
        vog_to_consider = vog_clusters_df[vog_clusters_df['vog_cluster'].isin(vog_cluster)]['vog']
        # st.dataframe(vog_to_consider)
        st.dataframe(VOG_GeneProduct.loc[VOG_GeneProduct['VOG'].isin(vog_to_consider.values), :].reset_index(drop=True))

# # st.subheader("Samples composition")
# # st.caption("Percentage of samples containing more than a specific percentage (50-60-70-80-90) of protein count")

# # but1, but2 = st.columns(2)
# # with but1:
# #        hide = st.button("Considered specie only")
# # with but2:
# #        show = st.button("All species")

# # # Function to display DataFrame based on button click
# # def display_data(df):
# #     if show == True:
# #         st.dataframe(df.set_index('species'))
# #     elif hide == True:
# #         st.dataframe(df[df['species']==option].set_index('species'))
# #     else:
# #         st.dataframe(df.set_index('species'))

# # # Call function to display data
# # display_data(mgss_candidates)
