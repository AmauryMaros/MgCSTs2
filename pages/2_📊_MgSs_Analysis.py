import streamlit as st
import pandas as pd

st.set_page_config(layout="centered")

# Genes based data
clust = pd.read_csv("Data/mgss.clustering.parameters.csv").rename(columns={"Unnamed: 0" : "species"})
mgss_coverage = pd.read_csv("Data/mgSs.coverage.stats.csv").rename(columns={"Unnamed: 0" : "sub_species", "as.factor(sample_cluster)" : "sample_cluster"})
species = mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]).unique()
mgss_candidates = pd.read_csv("Data/potential.mgss.candidates.csv")

# VOG based data
vog_clust = pd.read_csv("Data/mgss.clustering.parameters.vog.csv").rename(columns={"Unnamed: 0" : "species"})
vog_mgss_coverage = pd.read_csv("Data/vog.mgSs.coverage.stats.csv").rename(columns={"Unnamed: 0" : "sub_species", "as.factor(sample_cluster)" : "sample_cluster"})
vog_species = vog_mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]).unique()


st.title("Metagenomic Subspecies")

st.subheader("Visualizations")


a = len(species)
st.write(a)
option = st.selectbox("Species", species)
tab1, tab2, tab3 = st.tabs(["Species coverage", "Subspecies stats", "Presence Absence heatmap"])

with tab1:
        st.subheader("Genes based")
        col1, col2 = st.columns(2)
        with col1 :
        #     fig1 = Image.open("Medias/mgss_coverage_png/" + option + "_subspecies_coverage_boxplot.png") 
        #     st.image(fig1)
                st.image("Medias/mgss_coverage_png/" + option + "_subspecies_coverage_boxplot.png")

        with col2 :
        #     fig2 = Image.open("Medias/mgss_coverage_png/" + option + "_subspecies_coverage_by_NoGenes.png")
        #     st.image(fig2)
                st.image("Medias/mgss_coverage_png/" + option + "_subspecies_coverage_by_NoGenes.png")
        
        st.subheader("VOG based")
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
        st.subheader("Genes based")
        st.write("Clustering parameters")
        df2 = clust[clust['species'] == option].drop('Number_of_unassigned', axis = 1)
        st.dataframe(df2)

        st.write("Subspecies stats")
        st.dataframe(mgss_coverage[mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]) == option])

        st.subheader("VOG based")
        st.write("Clustering parameters")
        vog_df2 = vog_clust[vog_clust['species'] == option].drop('Number_of_unassigned', axis = 1)
        st.dataframe(vog_df2)

        st.write("Subspecies stats")
        st.dataframe(vog_mgss_coverage[vog_mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]) == option])



with tab3 :
        # image = Image.open("Medias/heatmap_presence_absence/_" + option + "_heatmap_presence_absence.png")
        # st.image(image)
        col1, col2 = st.columns(2)
        with col1 :
            st.subheader("Genes based")
            st.image("Medias/heatmap_presence_absence/_" + option + "_heatmap_presence_absence.png")
        with col2 :
            st.subheader("VOG based")
            st.image("Medias/vog_heatmap_presence_absence/_" + option + "_heatmap_presence_absence.png")


# st.subheader("Samples composition")
# st.caption("Percentage of samples containing more than a specific percentage (50-60-70-80-90) of protein count")

# but1, but2 = st.columns(2)
# with but1:
#        hide = st.button("Considered specie only")
# with but2:
#        show = st.button("All species")

# # Function to display DataFrame based on button click
# def display_data(df):
#     if show == True:
#         st.dataframe(df.set_index('species'))
#     elif hide == True:
#         st.dataframe(df[df['species']==option].set_index('species'))
#     else:
#         st.dataframe(df.set_index('species'))

# # Call function to display data
# display_data(mgss_candidates)
