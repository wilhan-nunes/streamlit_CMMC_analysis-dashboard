import streamlit as st


def render_welcome_message():
    st.markdown("""
     ## 🔬 Welcome to the CMMC Analysis Dashboard

     This interactive dashboard enables comprehensive analysis of **Collaborative Microbial Metabolite Center** enrichment workflow results combined with **Feature-Based Molecular Networking (FBMN)** data.
     [Access the full documentation](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/metaboapp_CMMC_dashboard/) for more details.
    """)

    st.info("""
    - This application is part of the GNPS downstream analysis ecosystem known as **MetaboApps**.
    - If you encounter any issues or have suggestions, please reach out to the app maintainers.
    - [Checkout other tools](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/metaboapps_overview/)
    """)

    st.markdown("""
     ### 📖 Citation
     Mannochio-Russo H, Nunes WDG et al. **Community Curation of Microbial Metabolites for The Mechanistic Analysis of Metabolomics Data.** Manuscript under preparation
     
     ### 🧭 **What This Tool Does**

     **Integrate Multiple GNPS2 workflow results:**
     - Combines CMMC enrichment results with FBMN quantification data
     - Merges metabolomics data with user-provided [sample metadata](#expected-metadata-table-format)
     - Creates merged tables and provides tools for comprehensive microbial metabolites analysis

     **Generate Interactive Visualizations:**
     - **Box Plots**: Compare metabolite abundances across different sample groups
     - **UpSet Plots**: Visualize overlaps and intersections in your data by metabolites source or origin
     - **Dynamic Filtering**: Explore data with flexible grouping and filtering options
     - **Advanced Network Visualization**: Visualize molecular networks with interactive features
     - **MicrobeMASST Search**: Search your spectra against a reference database of MS/MS data acquired from bacterial and fungal monocultures.

     ### 📘 **Getting Started**

     1. **Enter Task IDs**: Input your CMMC Enrichment and FBMN Task identifiers in the sidebar
     2. **Upload Metadata**: Provide your sample metadata table (CSV, Excel, TSV or TXT format)
     3. **Run Analysis**: Click the "🚀 Run Analysis" button to fetch and process your data
     4. **Explore Results**: Use the interactive plots to investigate your metabolomics data

     ### 🧩 **Analysis Features**
  
     **Box Plot Analysis:**
     - Compare metabolite abundances between sample groups
     - Filter by metadata attributes
     - Select specific features for detailed examination
     - Visualize statistical distributions and outliers

     **Metabolites Co-occurrence Visualization:**
     - Understand data intersections and unique elements
     - Group by source or origin classifications
     - Identify patterns in metabolite presence across samples
     
     **Molecular Network** (Advanced Visualization)
     - Visualize chemical relationships between metabolites based on spectral similarity
     - Explore networks of structurally related compounds across samples
     - Identify potential novel compounds or analogs through network propagation
     
     **MicrobeMASST Taxonomic Tree** (Advanced Visualization)
     - Map microbial metabolite origins onto a taxonomic hierarchy
     - Explore metabolite detection at different taxonomic levels (e.g., genus, species)

     ### **Tips for Best Results**

     - [Ensure your metadata table contains sample identifiers that match your FBMN data](#expected-metadata-table-format)
     - Use descriptive column names in your metadata for clearer visualizations
     - Try different grouping combinations to uncover hidden patterns in your data

     ---

     ### Ready to analyze your data?
     Configure your analysis parameters in the sidebar and click "🚀 Run Analysis" to begin!
     """)

    with st.expander("Sample Data Format"):
        st.subheader("Expected Metadata Table Format:")
        st.markdown("""
            | filename        | ATTRIBUTE_Treatment | ATTRIBUTE_Timepoint | ATTRIBUTE_Source | ATTRIBUTE_Origin     |
            |-----------------|---------------------|---------------------|------------------|----------------------|
            | Sample_01.mzML  | Control             | 24h                 | Plant            | Wild                 |
            | Sample_02.mzXML | Treated             | 24h                 | Plant            | Cultivated           |
            | Sample_03.mzML  | Control             | 48h                 | Fungal           | Wild                 |

         - **filename** (mandatory): Must match identifiers in your FBMN data
         - **Additional columns**: Can be used for grouping and filtering in analyses, and must start with "ATTRIBUTE_"
         - **File formats**: Accepts CSV, Excel, TSV, or TXT files
         """)
    with st.expander("🔗 About Task IDs"):
        st.markdown("""
         **CMMC Task ID**: Identifier from your GNPS2 CMMC enrichment analysis

         **FBMN Task ID**: Identifier from your GNPS2 Feature-Based Molecular Networking analysis

         Workflows to run those analysis are available at http://www.gnps2.org

         [FMBN Documentation](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/fbmn/) | 
         [CMMC Enrichment Documentation](https://cmmc.gnps2.org/network_enrichment/)

         """)
