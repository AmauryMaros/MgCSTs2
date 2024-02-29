import streamlit as st
import pandas as pd

st.set_page_config(layout="centered")

vog_clust = pd.read_csv("Data/mgss.clustering.parameters.vog.csv").rename(columns={"Unnamed: 0" : "species"})
vog_mgss_coverage = pd.read_csv("Data/vog.mgSs.coverage.stats.csv").rename(columns={"Unnamed: 0" : "sub_species", "as.factor(sample_cluster)" : "sample_cluster"})
vog_species = vog_mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]).unique()
# mgss_candidates = pd.read_csv("Data/potential.mgss.candidates.csv")


st.title("Metagenomic Subspecies")

st.subheader("Visualizations")


a = len(vog_species)
st.write(a)
option = st.selectbox("Species", vog_species)
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
        # image = Image.open("Medias/heatmap_presence_absence/_" + option + "_heatmap_presence_absence.png")
        # st.image(image)
        st.image("Medias/vog_heatmap_presence_absence/_" + option + "_heatmap_presence_absence.png")

