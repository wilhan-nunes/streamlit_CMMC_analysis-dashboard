import os

import streamlit as st
import streamlit.components.v1 as components

from microbemasst_search import run_microbemasst_search


def render_microbemasst_frame():
    st.subheader(":material/microbiology: MicrobeMASST Search")
    fbmn_task_id = st.session_state.get('fbmn_task_id', None)

    col1, col2, col3 = st.columns(3)

    with col1:
        usi_or_fid = st.text_input("Enter USI or Library ID:")
        prec_tol = st.number_input("Precursor m/z tolerance (ppm):", value=0.05, format="%.2f")
        mz_tol = st.number_input("m/z fragment tolerance (ppm):", value=0.05, format="%.2f")

    with col2:
        cos = st.number_input("Minimum cosine similarity:", value=0.7, min_value=0.0, max_value=1.0, format="%.2f")
        min_match_peaks = st.number_input("Minimum matched peaks:", value=3, min_value=1, step=1)

    with col3:
        analog_mass_below = st.number_input("Analog mass below (Da):", value=130, min_value=1, step=1)
        analog_mass_above = st.number_input("Analog mass above (Da):", value=140, min_value=1, step=1)
        use_analog = st.checkbox("Use Analog Masses", value=False)

    if st.button("Run Search"):
        if usi_or_fid.strip().lower().startswith("mzspec"):
            usi = usi_or_fid.strip()
        elif usi_or_fid.strip().isdigit():
            usi = f'mzspec:GNPS2:TASK-{fbmn_task_id}-nf_output/clustering/spectra_reformatted.mgf:scan:{usi_or_fid.strip()}'
        else:
            st.error(f"Input ':blue-badge[{usi_or_fid}]' is invalid. \n Please enter a valid USI or Library ID.")
            return

        with st.spinner("Running MicrobeMASST search..."):
            # Call the search function with the provided parameters
            out_path = run_microbemasst_search(
                usi, prec_tol, mz_tol, cos, min_match_peaks,
                analog_mass_below, analog_mass_above, use_analog
            )
        st.success(f"Search completed. Results saved in: {out_path}")
        all_files = os.listdir(out_path)
        print(all_files)

        if 'fastMASST_microbe.html' in all_files:

            with st.container(border=True):
                html_file = out_path + '/fastMASST_microbe.html'
                st.components.v1.html(
                    open(html_file).read(),
                    height=500,
                    scrolling=True)
                with open(html_file, 'rb') as f:
                    st.download_button(
                        label="Download tree HTML file",
                        data=f,
                        file_name='fastMASST_microbe.html',
                        mime='text/html',
                        type='tertiary',
                        icon=':material/download:'
                    )
        else:
            st.error("Error: The search didn't return any results.")
