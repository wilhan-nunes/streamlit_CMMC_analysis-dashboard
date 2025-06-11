import streamlit as st

import box_plot
import upset_plot
from utils import *


def main():
    st.set_page_config(
        page_title="CMMC Analysis Dashboard",
        page_icon="🦠",
        layout="wide"
    )

    # --- Query params for task IDs ---
    query_params = st.query_params
    default_cmmc_task_id = query_params.get("cmmc_task_id", "")
    default_fbmn_task_id = query_params.get("fbmn_task_id", "")

    st.title("🦠 CMMC Analysis Dashboard")
    st.markdown("---")

    # Sidebar configuration
    with st.sidebar:
        st.header("📊 Analysis Configuration")
        cmmc_task_id = st.text_input(
            "CMMC Enrichment Task ID",
            value=default_cmmc_task_id,
            placeholder="Enter CMMC Enrichment Task ID",
            help="Input your CMMC enrichment task identifier"
        )
        fbmn_task_id = st.text_input(
            "FBMN Task ID",
            value=default_fbmn_task_id,
            placeholder="Enter FBMN Task ID",
            help="Input your Feature-Based Molecular Network task identifier"
        )

        uploaded_file = st.file_uploader(
            "Upload Metadata Table",
            type=['csv', 'xlsx', 'tsv'],
            help="Upload your metadata table (CSV, Excel, or TSV format)"
        )

        # Display upload status
        if uploaded_file is not None:
            st.success(f"✅ File uploaded: {uploaded_file.name}")
            try:
                # Read the uploaded file
                if uploaded_file.name.endswith('.csv'):
                    loaded_metadata_df = pd.read_csv(uploaded_file)
                elif uploaded_file.name.endswith('.tsv'):
                    loaded_metadata_df = pd.read_csv(uploaded_file, sep='\t')
                else:  # Excel files
                    loaded_metadata_df = pd.read_excel(uploaded_file)

                st.session_state['metadata_df'] = loaded_metadata_df
                st.info(f"📋 Rows: {len(loaded_metadata_df)} | Columns: {len(loaded_metadata_df.columns)}")

                # Show preview
                with st.expander("Preview Data"):
                    st.dataframe(loaded_metadata_df.head(), use_container_width=True)

            except Exception as e:
                st.error(f"❌ Error reading file: {str(e)}")
        else:
            st.info("📤 Please upload a metadata table")

        st.markdown("---")

        # Analysis button
        run_analysis = st.button(
            "🚀 Run Analysis",
            type="primary",
            use_container_width=True,
            disabled=not (cmmc_task_id and fbmn_task_id and uploaded_file)
        )

        if st.button("Reset Analysis", type="secondary", use_container_width=True):
            st.session_state.clear()


    if run_analysis:
        st.session_state['run_analysis'] = True
        # Fetch enriched results and store in session state
        enriched_result = fetch_enriched_results(cmmc_task_id)
        st.session_state['enriched_result'] = enriched_result

        # fetch quantification data
        quant_file = fetch_file(fbmn_task_id, "quant_table.csv", "quant_table")
        if quant_file:
            df_quant = pd.read_csv(quant_file)
            st.session_state['df_quant'] = df_quant

    # Main content area
    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader("📈 UpSet Plot")
        # if run_analysis = true in session state, then generate upset plot

        if st.session_state.get('run_analysis'):
            group_by = st.segmented_control("Group by", ['Origin', "Source"], default="Origin")

            ss_enriched_result = st.session_state.get('enriched_result')
            upset_fig_source = upset_plot.generate_upset_plot(ss_enriched_result, by="source")
            upset_fig_origin = upset_plot.generate_upset_plot(ss_enriched_result, by="origin")

            if group_by == "Source":
                upset_fig = upset_fig_source
            else:
                upset_fig = upset_fig_origin

            st.pyplot(upset_fig, use_container_width=True)

        else:
            st.info("📋 Configure parameters in sidebar and run analysis to display UpSet plot")

    with col2:
        st.subheader("📊 Box Plots")
        if st.session_state.get('run_analysis'):
            groups = st.text_input('Input groups separated by commas', value="Stomach, Duodenum, Jejunum, Ileum, Cecum, Colon, Stool")
            formatted_groups = [i.strip() for i in groups.split(',')]

            quant = st.session_state.get('df_quant')
            metadata = st.session_state.get('metadata_df')
            ss_enriched_result = st.session_state.get('enriched_result')
            merged_data = box_plot.prepare_lcms_data(quant, metadata, ss_enriched_result)
            # from merged data create a input widget to select featureID (with input_namr) from merged_data.columns
            feat_id_dict = merged_data[['featureID', 'input_name']].drop_duplicates("featureID").set_index('featureID').to_dict(orient='index')
            feature_id = st.selectbox("Select Feature ID", [f"{k}: {v.get('input_name')}" for k,v in feat_id_dict.items()])

            # column_to_plot = st.selectbox("Column to plot", [val for val in merged_data.columns.tolist() if "ATTRIBUTE" in val], )

            st.plotly_chart(
                box_plot.plot_boxplots_by_group(
                    merged_data,
                    formatted_groups,
                    int(feature_id.split(":")[0]),
                    "ATTRIBUTE_UBERONBodyPartName"
                ),
                use_container_width=True
            )


    # Status information
    st.markdown("---")

    # Display current configuration
    with st.expander("🔧 Current Configuration", expanded=False):
        col1, col2, col3 = st.columns(3)

        with col1:
            st.write("**CMMC Task ID:**")
            st.code(cmmc_task_id if cmmc_task_id else "Not set")

        with col2:
            st.write("**FBMN Task ID:**")
            st.code(fbmn_task_id if fbmn_task_id else "Not set")

        with col3:
            st.write("**Metadata File:**")
            st.code(uploaded_file.name if uploaded_file else "Not uploaded")

if __name__ == "__main__":
    main()