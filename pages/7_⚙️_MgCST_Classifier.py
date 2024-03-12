import streamlit as st
import pandas as pd
from io import StringIO
import subprocess
import base64
import time
from time import sleep
from stqdm import stqdm

st.title("In progress ...")

st.subheader("Inputs", divider='gray')

summary_abundance = st.file_uploader("Import your summary.Abundance.txt file here")
if summary_abundance is not None:
    summary_abundance_df = pd.read_csv(summary_abundance)
    st.dataframe(summary_abundance_df.head())

summary_abundance_NR = st.file_uploader("Import your summary.NR.abundance.txt file here")
if summary_abundance_NR is not None:
    summary_abundance_NR_df = pd.read_csv(summary_abundance_NR)


st.subheader("Run classifier", divider='gray')

run = st.button("Run")

if run :

    # for i in stqdm(range(50), backend=False, frontend=True):
    #     sleep(0.5)

    st.subheader("Outputs", divider='gray')

    with st.spinner('The classifier is running...'):
        process = subprocess.Popen(["Rscript", "mgCSTs_heatmap.R"])

        for _ in stqdm(range(90)):
            sleep(1)

        result = process.communicate()
        st.success('Done!')
        st.toast('Classifier has successfully run', icon='✅')

        with open("Medias/mgCST_heatmap.pdf", "rb") as f:
                base64_pdf = base64.b64encode(f.read()).decode('utf-8')

                btn = st.download_button(
                        label="Download file",
                        data=f,
                        file_name="mgCST_heatmap.pdf",
                        mime="image/png"
                    )
                


# def run_r_script():
#     r_script_path = "mgCSTs_heatmap.R"  # Replace with the actual path to your R script
#     result = subprocess.Popen(["Rscript", r_script_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
#     while result.poll() is None:
#         time.sleep(0.1)
#         st.text(result.stdout.readline())

#     return result.stdout.read()

# def main():
#     st.title("Streamlit App with Live R Script Output")

#     if st.button("Run R Script"):
#         st.text("Running R script...")

#         with st.spinner("Please wait... Running R script."):
#             console_output = run_r_script()

#         st.text("R Script Output:")
#         st.text(console_output)

#         with open("Medias/mgCST_heatmap.pdf", "rb") as f:
#             base64_pdf = base64.b64encode(f.read()).decode('utf-8')
#             btn = st.download_button(
#             label="Download file",
#             data=f,
#             file_name="mgCST_heatmap.pdf",
#             mime="image/png"
#                     )

# if __name__ == "__main__":
#     main()
