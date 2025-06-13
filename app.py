import streamlit as st

import box_plot
import upset_plot
from utils import *


def main():
    st.set_page_config(
        page_title="CMMC Analysis Dashboard", page_icon="favicon.png", layout="centered"
    )

    # --- Query params for task IDs ---
    # http://localhost:8501/?cmmc_task_id=21c17a8de65041369d607493140a367f&fbmn_task_id=58e0e2959ec748049cb2c5f8bb8b87dc
    query_params = st.query_params
    default_cmmc_task_id = query_params.get("cmmc_task_id", "")
    default_fbmn_task_id = query_params.get("fbmn_task_id", "")

    # Sidebar configuration
    with st.sidebar:
        st.header("📊 Analysis Configuration")
        cmmc_task_id = st.text_input(
            "CMMC Enrichment Task ID",
            value=default_cmmc_task_id,
            placeholder="Enter CMMC Enrichment Task ID",
            help="Input your CMMC enrichment task identifier",
        )
        fbmn_task_id = st.text_input(
            "FBMN Task ID",
            value=default_fbmn_task_id,
            placeholder="Enter FBMN Task ID",
            help="Input your Feature-Based Molecular Network task identifier",
        )

        uploaded_file = st.file_uploader(
            "Upload Metadata Table",
            type=["csv", "xlsx", "tsv"],
            help="Upload your metadata table (CSV, Excel, or TSV format)",
        )

        # Display upload status
        if uploaded_file is not None:
            st.success(f"✅ File uploaded: {uploaded_file.name}")
            try:
                # Read the uploaded file
                if uploaded_file.name.endswith(".csv"):
                    loaded_metadata_df = pd.read_csv(uploaded_file)
                elif uploaded_file.name.endswith(".tsv"):
                    loaded_metadata_df = pd.read_csv(uploaded_file, sep="\t")
                else:  # Excel files
                    loaded_metadata_df = pd.read_excel(uploaded_file)

                st.session_state["metadata_df"] = loaded_metadata_df
                st.info(
                    f"📋 Rows: {len(loaded_metadata_df)} | Columns: {len(loaded_metadata_df.columns)}"
                )

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
            disabled=not (cmmc_task_id and fbmn_task_id and uploaded_file),
        )

        if st.button("Reset Analysis", type="secondary", use_container_width=True):
            st.session_state.clear()
            st.rerun()

    if run_analysis:
        st.session_state["run_analysis"] = True
        # Fetch enriched results and store in session state
        enriched_result = fetch_enriched_results(cmmc_task_id)
        enriched_result['input_molecule_origin'] = enriched_result['input_molecule_origin'].str.replace(
            ' (e.g., natural products and other specialized metabolites)', '')

        st.session_state["enriched_result"] = enriched_result

        # fetch quantification data
        quant_file = fetch_file(fbmn_task_id, "quant_table.csv", "quant_table")
        if quant_file:
            df_quant = pd.read_csv(quant_file)
            st.session_state["df_quant"] = df_quant

        st.session_state["merged_df"] = box_plot.prepare_lcms_data(
            df_quant, loaded_metadata_df, enriched_result
        )

    # Initial page loaded if "run_analysis" not in st.session_state
    if not st.session_state.get('run_analysis'):
        # Welcome page content
        from welcome import render_welcome_message
        render_welcome_message()

    # Main content area
    if st.session_state.get("run_analysis"):

        st.title("🦠 CMMC Analysis Dashboard")
        st.markdown("---")

        st.subheader("Data Overview",
                     help="Select the column that contains the groups you want to compare to see a boxplot for each detected feature.")
        data_overview_df = st.session_state.get("merged_df")

        col1, col2 = st.columns(2)
        with col1:
            column_select = st.selectbox(
                "Select column", [i for i in data_overview_df.columns]
            )
            if column_select:
                group_by = st.multiselect(
                    "Select groups to compare",
                    [i for i in data_overview_df[column_select].unique()],
                    key="a",
                )
        with col2:
            # Filter data_overview_df based on the selected column and value
            input1, input2 = st.columns(2)
            with input1:
                first = st.selectbox('Column', ["input_molecule_origin", "input_source"])
            with input2:
                origin_list = ['Ambiguous',
                               'De novo biosynthesis by microbes',
                               'Diet', 'Drug', 'Exposure', 'Exposure/diet', 'Host',
                               'Host metabolism of microbial metabolites', 'Insecticides/pesticides',
                               'Microbial metabolism of drugs', 'Microbial metabolism of food molecules',
                               'Microbial metabolism of host-derived molecules',
                               'Microbial metabolism of microbial-derived molecules',
                               'Microbial metabolism of other human-made molecules', 'Unknown/Undefined']
                source_list = ['Microbial', 'Host', 'Diet', 'Unknown', 'Ambiguous', 'Drug', 'Exposure',
                               'Pesticides/insecticides', 'Other human-made molecules']
                second = st.multiselect("Value", origin_list if first == "input_molecule_origin" else source_list)

            data_overview_df['input_clean'] = (
                data_overview_df[first]
                .fillna('')
                .str.replace(r'\s+and\s+', ';', regex=True)
                .str.split(';')
                .apply(lambda items: list({item.strip().lower() for item in items if item}))
            )

            if st.checkbox("Use column and value filters"):
                target_set = set([i.lower() for i in second])
                mask = data_overview_df['input_clean'].apply(
                    lambda x: bool(set(x) & target_set) if isinstance(x, list) else False)
                data_overview_df = data_overview_df[mask]
            # data_overview_df = data_overview_df[data_overview_df['input_clean'] == set(second)]

            feat_id_dict = (
                data_overview_df[["featureID", "input_name"]]
                .drop_duplicates("featureID")
                .set_index("featureID")
                .to_dict(orient="index")
            )

            fid_items = [f"{k}: {v.get('input_name')}" for k, v in feat_id_dict.items()]
            feature_id = st.selectbox(
                f"Select Feature ID :blue-badge[{len(fid_items)} item(s)]",
                fid_items,
                key="b",
            )
        st.plotly_chart(
            box_plot.plot_boxplots_by_group(
                data_overview_df,
                groups1=group_by,  # this will be on x axis
                column1=column_select,
                feature_id=int(feature_id.split(":")[0]),
            ),
            use_container_width=True,
            key="graph1",
        )

    if st.session_state.get("run_analysis"):
        st.subheader("📊 Box Plots")
        quant = st.session_state.get("df_quant")
        metadata = st.session_state.get("metadata_df")
        ss_enriched_result = st.session_state.get("enriched_result")
        merged_data = box_plot.prepare_lcms_data(
            quant, metadata, ss_enriched_result
        )

        col_attr1, col_attr2 = st.columns(2)
        with col_attr1:
            selected_attribute1 = st.selectbox(
                "Metadata prefilter (group 1)",
                [i for i in metadata.columns],
                help="Prefilter the data based on the given group (optional)",
            )
        with col_attr2:
            selected_attribute2 = st.selectbox(
                "Metadata group 2 (x axis)", [i for i in metadata.columns]
            )

        col1, col2 = st.columns([1, 1])
        with col1:
            if selected_attribute1:
                groups1 = st.selectbox(
                    "Group 1 (single)",
                    [i for i in metadata[selected_attribute1].unique()],
                )
            else:
                groups1 = None

        with col2:
            if selected_attribute2:
                groups2 = st.multiselect(
                    "Group 2 (multi)",
                    [i for i in metadata[selected_attribute2].unique()],
                    accept_new_options=True,
                )

        # from merged data create an input widget to select featureID (with input_name) from merged_data.columns
        feat_id_dict = (
            merged_data[["featureID", "input_name"]]
            .drop_duplicates("featureID")
            .set_index("featureID")
            .to_dict(orient="index")
        )
        feature_id = st.selectbox(
            "Select Feature ID",
            [f"{k}: {v.get('input_name')}" for k, v in feat_id_dict.items()],
        )

        prefilter = selected_attribute1 if selected_attribute1 != "None" else None
        try:
            st.plotly_chart(
                box_plot.plot_boxplots_by_group(
                    merged_data,
                    groups2,
                    [groups1],
                    int(feature_id.split(":")[0]),
                    selected_attribute2,
                    prefilter,
                ),
                use_container_width=True,
                key="graph2",
            )
        except:
            st.warning("Select all required fields to see the boxplot")

    if st.session_state.get("run_analysis"):
        st.subheader("📈 UpSet Plot")
        group_by = st.segmented_control(
            "Group by", ["Source", "Origin"], default="Source"
        )

        ss_enriched_result = st.session_state.get("enriched_result")
        upset_fig_source = upset_plot.generate_upset_plot(
            ss_enriched_result, by="source"
        )
        upset_fig_origin = upset_plot.generate_upset_plot(
            ss_enriched_result, by="origin"
        )

        if group_by == "Source":
            upset_fig = upset_fig_source
        else:
            upset_fig = upset_fig_origin

        st.pyplot(upset_fig, use_container_width=False)


if __name__ == "__main__":
    main()
