import os

import streamlit as st
import streamlit.components.v1 as components

from microbemasst_search import run_microbemasst_search


def render_microbemasst_frame():
    st.subheader("MicrobeMASST Search")

    col1, col2, col3 = st.columns(3)
    with col1:
        usi = st.text_input("Enter USI or Library ID:")
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
                st.components.v1.html(
                    open(out_path + '/fastMASST_microbe.html').read(),
                    height=800,
                    scrolling=True)
        else:
            st.error("Error: The search didn't return any results.")
