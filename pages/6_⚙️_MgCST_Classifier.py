import streamlit as st
import pandas as pd
import base64

st.set_page_config(layout="wide")

######################################## ----- DATA IMPORTATION ----- ########################################

relabund = pd.read_csv("Data/vog_relabund_w_mgCSTs_12Mar2024.csv")
mgcsts = pd.read_csv("Data/vog_mgCSTs_12Mar2024.csv")
vog_abundance = pd.read_csv("Data/vog_norm_counts_mgSs_mgCST_12Mar2024.csv")

######################################## ----- DATA PROCESSING ----- ########################################

relabund = relabund.rename(columns = {"Unnamed: 0":"sampleID"})
mgcsts = mgcsts.rename(columns = {"Unnamed: 0":"sampleID"})
vog_abundance = vog_abundance.rename(columns = {"Unnamed: 0":"sampleID"})


######################################## ----- PAGE CONTENT ----- ########################################

st.header("MgCSTs Classifier Outputs", divider = 'grey')

def displayPDF(file):
    # Opening file from file path
    with open(file, "rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')
    # Embedding PDF in HTML
    pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="300" height="400" type="application/pdf"></iframe>'
    # Displaying File
    st.markdown(pdf_display, unsafe_allow_html=True)

st.container()
col1 ,col2 = st.columns(2)
with col1 :
    st.subheader("Heatmap")
    displayPDF("Medias/vog_mgCST_heatmap_13Mar2024.pdf")

with col2 :
    st.subheader("Relative abundances of metagenomic subspecies (columns) by samples (rows)")
    st.dataframe(relabund)

st.container()
col3 ,col4 = st.columns(2)
with col3 :
    st.subheader("SampleID - MgCST assignment")
    st.dataframe(mgcsts)
with col4 :
    st.subheader("VOG abundance")
    st.dataframe(vog_abundance)
