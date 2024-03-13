import streamlit as st
import base64

def displayPDF(file):
    # Opening file from file path
    with open(file, "rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')

    # Embedding PDF in HTML
    pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="300" height="400" type="application/pdf"></iframe>'

    # Displaying File
    st.markdown(pdf_display, unsafe_allow_html=True)

col1, col2 = st.columns(2)
with col1 : 
    st.subheader("Gene based")
    displayPDF("Medias/mgCST_heatmap_13Mar2024.pdf")
with col2 :
    st.subheader("VOG based")
    displayPDF("Medias/vog_mgCST_heatmap_13Mar2024.pdf")