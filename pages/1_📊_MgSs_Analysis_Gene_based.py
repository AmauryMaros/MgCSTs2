import streamlit as st
import pandas as pd

st.set_page_config(layout="centered")

######################################## ----- DATA IMPORTATION ----- ########################################

clust = pd.read_csv("Data/mgss.clustering.parameters.csv").rename(columns={"Unnamed: 0" : "species"})
mgss_coverage = pd.read_csv("Data/mgSs.coverage.stats.csv").rename(columns={"Unnamed: 0" : "sub_species", "as.factor(sample_cluster)" : "sample_cluster"})


######################################## ----- DATA PROCESSING ----- ########################################

species = mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]).unique()
not_to_cluster = ["Alterileibacterium", "Anaerococcus", "Bacteroides", "Campylobacter", 
                  "Corynebacterium", "Gardnerella", "Gulosibacter", "Lactobacillus", 
                  "Limosilactobacillus", "MultiGenera", "Porphyromonas", "Prevotella", "Streptococcus"]

species = [s for s in species if s not in not_to_cluster]


######################################## ----- PAGE CONTENT ----- ########################################

st.title("Metagenomic Subspecies")

st.subheader("Visualizations")

option = st.selectbox("Species", species)

tab1, tab2, tab3 = st.tabs(["Species coverage", "Subspecies stats", "Presence Absence heatmap"])

with tab1:

        col1, col2 = st.columns(2)

        with col1 :
                st.image("Medias/mgss_coverage_png/" + option + "_subspecies_coverage_boxplot.png")

        with col2 :
                st.image("Medias/mgss_coverage_png/" + option + "_subspecies_coverage_by_NoGenes.png")

with tab2 :

        st.subheader("Clustering parameters")
        df2 = clust[clust['species'] == option].drop('Number_of_unassigned', axis = 1)
        st.dataframe(df2)

        st.subheader("Subspecies stats")
        st.dataframe(mgss_coverage[mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]) == option])

with tab3 :

        st.image("Medias/heatmap_presence_absence/_" + option + "_heatmap_presence_absence.png")

