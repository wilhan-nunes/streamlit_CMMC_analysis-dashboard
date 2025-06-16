import streamlit as st


def render_welcome_message():
    st.markdown("""
     ## 🔬 Welcome to the CMMC Analysis Dashboard

     This interactive dashboard enables comprehensive analysis of **Collaborative Microbial Metabolite Center** enrichment results combined with **Feature-Based Molecular Networking (FBMN)** data.

     ### 🎯 **What This Tool Does**

     **Integrate Multiple GNPS2 workflow results:**
     - Combines CMMC enrichment results with FBMN quantification data
     - Merges metabolomics data with user-provided [metadata](#expected-metadata-table-format)
     - Creates unified datasets for comprehensive microbial metabolites analysis

     **Generate Interactive Visualizations:**
     - **Box Plots**: Compare metabolite abundances across different sample groups
     - **UpSet Plots**: Visualize overlaps and intersections in your data by metabolites source or origin
     - **Dynamic Filtering**: Explore data with flexible grouping and filtering options

     ### **Getting Started**

     1. **Enter Task IDs**: Input your CMMC Enrichment and FBMN Task identifiers in the sidebar
     2. **Upload Metadata**: Provide your sample metadata table (CSV, Excel, or TSV format)
     3. **Run Analysis**: Click the "🚀 Run Analysis" button to fetch and process your data
     4. **Explore Results**: Use the interactive plots to investigate your metabolomics data

     ### **Analysis Features**

     **Overview Box Plot Analysis:**
     - Allows a broader overview of the data
     - Select each metadata attribute column and value to use as groups
     - Filter microbial metabolites by source or origin 
     - Select specific features for detailed examination
     - Download all features plots as a pdf file
     
     **Box Plot Analysis:**
     - Compare metabolite abundances between sample groups
     - Filter by metadata attributes
     - Select specific features for detailed examination
     - Visualize statistical distributions and outliers
     - Download all features plots as a pdf file

     **UpSet Plot Visualization:**
     - Understand data intersections and unique elements
     - Group by source or origin classifications
     - Identify patterns in metabolite presence across samples

     ### **Tips for Best Results**

     - [Ensure your metadata table contains sample identifiers that match your FBMN data](#expected-metadata-table-format)
     - Use descriptive column names in your metadata for clearer visualizations
     - Try different grouping combinations to uncover hidden patterns in your data

     ---

     **Ready to analyze your data?** Configure your analysis parameters in the sidebar and click "🚀 Run Analysis" to begin!
     """)
    # Optional: Add some sample data info or examples
    with st.expander("Sample Data Format"):
        st.subheader("Expected Metadata Table Format:")
        st.markdown("""
         | filename | Treatment | Timepoint | Source | Origin |
         |-----------|-----------|-----------|---------|---------|
         | Sample_01.mzML | Control   | 24h       | Plant   | Wild    |
         | Sample_02.mzXML | Treated   | 24h       | Plant   | Cultivated |
         | Sample_03.mzML | Control   | 48h       | Fungal  | Wild    |

         - **filename** (mandatory): Must match identifiers in your FBMN data
         - **Additional columns**: Can be used for grouping and filtering in analyses
         """)
    with st.expander("🔗 About Task IDs"):
        st.markdown("""
         **CMMC Task ID**: Identifier from your GNPS2 CMMC enrichment analysis

         **FBMN Task ID**: Identifier from your GNPS2 Feature-Based Molecular Networking analysis

         Workflows to run those analysis are available at http://www.gnps2.org

         [FMBN Documentation](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/fbmn/) | 
         [CMMC Enrichment Documentation](https://cmmc.gnps2.org/network_enrichment/)

         """)