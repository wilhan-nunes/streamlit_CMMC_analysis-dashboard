from streamlit.components.v1 import html

import upset_plot
from enhanced_boxplot import render_statistical_boxplot_tab
from network_cluster_plotter import *
from utils import *
from gnpsdata import workflow_fbmn


@st.cache_data
def get_gnps2_fbmn_metadata_table(taskid):
    return workflow_fbmn.get_metadata_dataframe(taskid, gnps2=True)


def render_sidebar():
    global cmmc_task_id, fbmn_task_id, uploaded_quant_file, loaded_metadata_df, include_all_features, run_analysis
    with st.sidebar:
        st.header("📊 Analysis Configuration")
        load_example = st.checkbox("Load Example Data", value=False, key="load_example")
        if not load_example:
            col_cmmc_input, col_cmmc_link = st.columns([8, 1])
            with col_cmmc_input:
                cmmc_task_id = st.text_input(
                    "CMMC Enrichment Task ID",
                    value=default_cmmc_task_id,
                    placeholder="Enter CMMC Enrichment Task ID",
                    help="Input your CMMC enrichment task identifier",
                    key="cmmc_task_id",
                )
            with col_cmmc_link:
                if cmmc_task_id and len(cmmc_task_id) == 32:
                    st.markdown('<div style="height: 28px;"></div>', unsafe_allow_html=True)
                    st.link_button("", f"https://gnps2.org/status?task={cmmc_task_id}", icon=':material/arrow_outward:', use_container_width=True)
            validate_task_id_input(cmmc_task_id, validation_str="cmmc")

            col_fbmn_input, col_fbmn_link = st.columns([8, 1])
            with col_fbmn_input:
                fbmn_task_id = st.text_input(
                    "FBMN Task ID",
                    value=default_fbmn_task_id,
                    placeholder="Enter FBMN Task ID",
                    help="Input your Feature-Based Molecular Network task identifier",
                    key="fbmn_task_id",
                )
            with col_fbmn_link:
                if fbmn_task_id and len(fbmn_task_id) == 32:
                    st.markdown('<div style="height: 28px;"></div>', unsafe_allow_html=True)
                    st.link_button("", f"https://gnps2.org/status?task={fbmn_task_id}", icon=':material/arrow_outward:',  use_container_width=True)
            validate_task_id_input(fbmn_task_id, 'feature_based')

            # Check if metadata is available from FBMN task
            fbmn_metadata_available = False
            if fbmn_task_id and len(fbmn_task_id) == 32:
                try:
                    metadata_from_fbmn = get_gnps2_fbmn_metadata_table(fbmn_task_id)
                    if isinstance(metadata_from_fbmn, pd.DataFrame) and not metadata_from_fbmn.empty:
                        fbmn_metadata_available = True
                        st.info("✓ Metadata table found in FBMN task", icon=":material/info:")
                except Exception:
                    fbmn_metadata_available = False

            # Metadata upload section
            if fbmn_metadata_available:
                use_fbmn_metadata = st.checkbox(
                    "Use metadata from FBMN task",
                    value=True,
                    help="Use the metadata table from your FBMN task, or upload a different one below",
                    key="use_fbmn_metadata"
                )
                
                if use_fbmn_metadata:
                    loaded_metadata_df = metadata_from_fbmn
                    st.session_state["metadata_df"] = loaded_metadata_df
                    st.success(
                        f"Rows: {len(loaded_metadata_df)} | Columns: {len(loaded_metadata_df.columns)}",
                        icon=":material/task:",
                    )
                    
                    # Show preview
                    with st.expander("Preview Data", icon=":material/visibility:"):
                        st.dataframe(loaded_metadata_df.head(), use_container_width=True)
                    
                    uploaded_metadata_file = True  # Set flag to enable analysis button
                else:
                    uploaded_metadata_file = st.file_uploader(
                        "Upload Metadata Table",
                        type=["csv", "xlsx", "tsv", "txt"],
                        help="Upload your metadata table (CSV, Excel, TSV or TXT format)",
                    )
            else:
                uploaded_metadata_file = st.file_uploader(
                    "Upload Metadata Table",
                    type=["csv", "xlsx", "tsv", "txt"],
                    help="Upload your metadata table (CSV, Excel, TSV or TXT format)",
                )

            if st.checkbox("Use uploaded quantification table", key="use_quant_table"):
                uploaded_quant_file = st.file_uploader(
                    "Upload Quantification Table",
                    type=["csv", "xlsx", "tsv", "txt"],
                    help="Upload your quantification table (CSV, Excel, TSV or TXT format)",
                )

            # Display upload status
            if uploaded_metadata_file is not None and uploaded_metadata_file is not True:
                try:
                    # Read the uploaded file
                    loaded_metadata_df = load_uploaded_file_df(uploaded_metadata_file)

                    if "filename" not in loaded_metadata_df.columns:
                        st.warning(
                            "Your metadata file must contain a 'filename' column",
                            icon=":material/warning:",
                        )
                    st.session_state["metadata_df"] = loaded_metadata_df
                    st.success(
                        f"Rows: {len(loaded_metadata_df)} | Columns: {len(loaded_metadata_df.columns)}",
                        icon=":material/task:",
                    )

                    metadata_cols = loaded_metadata_df.columns.tolist()
                    attrib_cols = [col for col in metadata_cols if "ATTRIBUTE_" in col]
                    if not attrib_cols:
                        st.warning(
                            "No columns with 'ATTRIBUTE_' prefix found. Select columns below to add the prefix.",
                            icon=":material/warning:",
                        )
                        selected_cols = st.multiselect(
                            "Select columns to add 'ATTRIBUTE_' prefix:",
                            [col for col in metadata_cols if not col.startswith("ATTRIBUTE_")],
                            help="Choose columns to be renamed with 'ATTRIBUTE_' prefix."
                        )
                        if selected_cols:
                            loaded_metadata_df.rename(
                                columns={col: f"ATTRIBUTE_{col}" for col in selected_cols},
                                inplace=True
                            )
                            st.success(
                                f"Added 'ATTRIBUTE_' prefix to: {', '.join(selected_cols)}",
                                icon=":material/task:",
                            )

                    # Show preview
                    with st.expander("Preview Data", icon=":material/visibility:"):
                        st.dataframe(loaded_metadata_df.head(), use_container_width=True)

                except Exception as e:
                    st.error(f"Error reading file: {str(e)}", icon=":material/error:")
            elif uploaded_metadata_file is None:
                st.info("📤 Please upload a metadata table")
        else:
            # loads Quinn's 2020 example data https://doi.org/10.1038/s41586-020-2047-9
            cmmc_task_id = "7f53b63490c945e980dfa10273a296cd"
            validate_task_id_input(cmmc_task_id, validation_str="cmmc")
            fbmn_task_id = "58e0e2959ec748049cb2c5f8bb8b87dc"
            st.session_state['fbmn_task_id'] = fbmn_task_id
            st.session_state['cmmc_task_id'] = cmmc_task_id
            uploaded_metadata_file = open('data/metadata_quinn2020.tsv', 'rb')

            st.write("Using example data from Quinn et al. 2020: https://doi.org/10.1038/s41586-020-2047-9")
            st.write(f"CMMC Task ID: [**{cmmc_task_id}**](https://gnps2.org/status?task={cmmc_task_id})")
            st.write(f"FBMN Task ID: [**{fbmn_task_id}**](https://gnps2.org/status?task={fbmn_task_id})")

            loaded_metadata_df = load_uploaded_file_df(uploaded_metadata_file)
            st.session_state["metadata_df"] = loaded_metadata_df
            with st.expander("Preview Data", icon=":material/visibility:"):
                st.dataframe(loaded_metadata_df.head(), use_container_width=True)

        st.markdown("---")
        include_all_features = st.checkbox(
            "Include all features in the analysis",
            value=False,
            help="If unchecked, only features with CMMC matches will be shown",
            key="include_all_features",
        )
        if include_all_features:
            st.warning(
                "This option will include all features in the analysis, even those without CMMC matches. "
                "This may lead to a larger dataset and longer processing time.",
                icon=":material/warning:",
            )
        # Analysis button
        run_analysis = st.button(
            "🚀 Run Analysis",
            type="primary",
            use_container_width=True,
            disabled=not (cmmc_task_id and fbmn_task_id and uploaded_metadata_file),
        )

        if st.button("Reset Analysis", type="secondary", use_container_width=True):
            st.session_state.clear()
            st.session_state["load_example"] = False
            st.rerun()

        st.subheader("Contributors")
        st.write("""<a href="https://sites.google.com/view/helenamrusso/home" target="_blank">Helena Russo PhD</a> - UC San Diego<br>
<a href="https://scholar.google.com/citations?user=4cPVoeIAAAAJ" target="_blank">Wilhan Nunes PhD</a> - UC San Diego
""", unsafe_allow_html=True)

        st.subheader("Documentations and Resources")
        st.write("""<a href="https://cmmc.gnps2.org/network_enrichment/">CMMC Enrichment Workflow</a><br>
                <a href="https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/fbmn/">Feature Based Molecular Networking</a><br>
                <a href="https://cmmc-kb.gnps2.org" target="_blank">CMMC Knowledge Base</a>""",
                 unsafe_allow_html=True)


# Set page configuration
# TODO: Bump version
app_version = "2025-09-24"
page_title = "CMMC Analysis Dashboard"
git_hash = get_git_short_rev()
repo_link = "https://github.com/wilhan-nunes/streamlit_CMMC_analysis-dashboard"

menu_items={"about": f"**CMMC Dashboard version: {app_version}**"}
st.set_page_config(
    page_title=page_title, page_icon="favicon.png", layout="wide", menu_items={"About": (f"**App version**: {app_version} | "
                          f"[**Git Hash**: {git_hash}]({repo_link}/commit/{git_hash})")}
)

html('<script async defer data-website-id="74bc9983-13c4-4da0-89ae-b78209c13aaf" src="https://analytics.gnps2.org/umami.js"></script>', width=0, height=0)
html('<script async defer data-website-id="1a70457c-e452-4e8f-822b-9983d742e174" src="https://analytics.gnps2.org/umami.js"></script>', width=0, height=0)

# --- Query params for task IDs ---
# http://localhost:8501/?cmmc_task_id=21c17a8de65041369d607493140a367f&fbmn_task_id=58e0e2959ec748049cb2c5f8bb8b87dc
query_params = st.query_params
default_cmmc_task_id = query_params.get("cmmc_task_id", "")
default_fbmn_task_id = query_params.get("fbmn_task_id", "")

# Sidebar configuration
render_sidebar()


def _process_data():
    print("Processing triggered...")
    # Create progress indicators
    progress_bar = st.progress(0)
    status_text = st.empty()
    try:
        # Step 1: Fetch enriched results
        status_text.text("Fetching CMMC enrichment results...")
        progress_bar.progress(10)

        enriched_result = fetch_enriched_results(cmmc_task_id)
        enriched_result["input_molecule_origin"] = enriched_result[
            "input_molecule_origin"
        ].apply(
            lambda x: str(x).replace(
                " (e.g., natural products and other specialized metabolites)", ""
            )
        )
        
        # Standardize source and origin columns for UpSet plot
        enriched_result["input_source_clean"] = (
            enriched_result["input_source"]
            .fillna("")
            .str.replace(r"\s+and\s+", ";", regex=True)
            .str.split(";")
            .apply(lambda items: list({item.strip() for item in items if item}))
        )
        enriched_result["input_molecule_origin_clean"] = (
            enriched_result["input_molecule_origin"]
            .fillna("")
            .str.replace(r"\s+and\s+", ";", regex=True)
            .str.split(";")
            .apply(lambda items: list({item.strip() for item in items if item}))
        )
        
        progress_bar.progress(25)

        # Step 2: Fetch quantification data
        status_text.text("Fetching quantification data...")
        include_all_features = st.session_state.get("include_all_features", False)

        if not st.session_state.get("use_quant_table", False):
            st.toast(
                "No quantification table uploaded. Using quantification table from FBMN job.",
                icon=":material/data_info_alert:",
            )

            quant_file = fbmn_quant_download_wrapper(fbmn_task_id)
        else:
            quant_file = load_uploaded_file_df(uploaded_quant_file)

        if quant_file is None:
            st.error("Failed to fetch quantification data")
            st.stop()

        progress_bar.progress(40)

        # Step 3: Process quantification data
        status_text.text("Processing quantification data...")
        if isinstance(quant_file, str):  # If it's a file path
            df_quant = pd.read_csv(quant_file)
        else:  # If it's already a DataFrame
            df_quant = quant_file


        progress_bar.progress(55)

        # Step 4: Merge data (this is the potentially slow step)
        if include_all_features:
            status_text.text("Merging all features (this may take a while for large datasets)...")
            st.info("Processing all features - this may take some time depending on dataset size.",
                    icon=":material/hourglass_top:")
        else:
            status_text.text("Merging CMMC-matched features...")

        progress_bar.progress(70)

        # Use the optimized function
        merged_df = prepare_lcms_data(
            df_quant, loaded_metadata_df, enriched_result, include_all_features
        )


        progress_bar.progress(85)

        # Step 5: Fetch network data
        status_text.text("Fetching molecular network data...")
        graphml_file_name = fetch_cmmc_graphml(
            cmmc_task_id
        )

        progress_bar.progress(90)

        # Step 6: Generate UpSet plots
        status_text.text("Generating UpSet plots...")
        try:
            upset_fig_source = upset_plot.generate_upset_plot(
                enriched_result, by="source"
            )
            upset_fig_origin = upset_plot.generate_upset_plot(
                enriched_result, by="origin"
            )
            st.session_state["upset_fig_source"] = upset_fig_source
            st.session_state["upset_fig_origin"] = upset_fig_origin
        except Exception as e:
            st.warning(f"Failed to generate UpSet plots: {str(e)}")
            st.session_state["upset_fig_source"] = None
            st.session_state["upset_fig_origin"] = None

        progress_bar.progress(100)

        # Success message
        status_text.text("✅ Analysis completed successfully!")

        # Clear progress indicators after a short delay
        import time
        time.sleep(1)
        progress_bar.empty()
        status_text.empty()

        st.session_state["run_analysis"] = True
        st.session_state["enriched_result"] = enriched_result
        st.session_state["df_quant"] = df_quant
        st.session_state["merged_df"] = merged_df
        st.session_state['G'] = nx.read_graphml(graphml_file_name)

    except Exception as e:
        progress_bar.empty()
        status_text.empty()
        st.error(f"Error during analysis: {str(e)}", icon=":material/error:")
        st.session_state["run_analysis"] = False
        st.stop()


if run_analysis:
    _process_data()

# Initial page loaded if "run_analysis" not in st.session_state
if not st.session_state.get("run_analysis"):
    # Welcome page content
    from welcome import render_welcome_message

    with st.container():
        render_welcome_message()

# Main content area
if st.session_state.get("run_analysis"):
    st.title("🦠 CMMC Analysis Dashboard")
    tabs = st.tabs(["🔭 Data Explorer", '⚙️ Advanced Visualizations'])
    with tabs[0]:
        # Box plot module
        merged_df = st.session_state.merged_df.infer_objects()
        render_statistical_boxplot_tab(merged_df, cmmc_task_id)

        st.markdown("---")

        # UpsetPlot module
        st.subheader(":green[:material/hub:] UpSet Plot")
        
        # Check if UpSet plots are available in session state
        upset_fig_source = st.session_state.get("upset_fig_source")
        upset_fig_origin = st.session_state.get("upset_fig_origin")
        
        if upset_fig_source is None and upset_fig_origin is None:
            st.error("UpSet plots are not available. Please re-run the analysis.")
        else:
            group_by = st.segmented_control(
                "Group metabolites by:", ["Source", "Origin"], default="Source"
            )

            if group_by == "Source":
                _, plot_col, _ = st.columns([1, 2, 1])
                upset_fig = upset_fig_source
                plot_type = "source"
            else:
                _, plot_col, _ = st.columns([1, 4, 1])
                upset_fig = upset_fig_origin
                plot_type = "origin"

            with plot_col:
                if upset_fig is not None:
                    st.image(upset_fig, use_container_width=False)
                    st.download_button(
                        label=":material/download: Download as SVG",
                        data=upset_fig,
                        file_name=f"upset_plot_{plot_type}.svg",
                        mime="image/svg+xml",
                        key='upset_plot_download'
                    )
                else:
                    st.error(f"UpSet plot for {plot_type} is not available. Please re-run the analysis.")

        # insert an expander card explaining how to interpret the upset plot
        with st.expander("How to interpret the UpSet plot", expanded=False):
            st.markdown(
                """
                The UpSet plot shows the co-occurrence of microbial metabolites across different sources or origins. 
                Each bar represents a unique combination of sources or origins, and the height of the bar indicates the number of features that match that combination.
                
                - **Left side**: The sets (sources or origins) that are being compared.
                - **Bottom side**: The intersections of these sets.
                - **Bars**: The height of each bar indicates how many features are present in that intersection.
                
                Use this plot to identify which sources or origins share the most microbial metabolites.
                """
            )
        st.markdown("---")


    with tabs[1]:
        # SETUP
        enriched_result = st.session_state.get("enriched_result")
        G = st.session_state['G']

        # Create a color_mapping from feature ID to component, filtering out single nodes in one pass
        nodes_dict = {
            str(row["query_scan"]): G.nodes[str(row["query_scan"])].get("component")
            for _, row in enriched_result.iterrows()
        }
        valid_nodes = {k: v for k, v in nodes_dict.items() if v != -1}

        # Build feature ID to name dict only for valid nodes
        feat_id_dict = {
            str(row["query_scan"]): row["input_name"]
            for _, row in enriched_result.iterrows()
            if str(row["query_scan"]) in valid_nodes
        }

        fid_labels = [
            f"{k}: {v} | Network {valid_nodes[k]}" for k, v in feat_id_dict.items()
        ]

        st.subheader("🕸️ Molecular Network Visualization")

        col_select, col_radio, col_deltas = st.columns([1, 1, 1])
        with col_select:
            selected_feature = st.selectbox(
                "Feature ID (no single nodes)",
                fid_labels,
                help="Annotated features that appear as single nodes in the network are excluded from this list.",
                width=500,
            )
        with col_radio:
            node_info = st.radio(
                "Node Legend", ["Feature ID", "Precursor m/z"], horizontal=True
            )
        with col_deltas:
            st.write('<div style="height: 40px;"></div>', unsafe_allow_html=True) # spacer
            show_deltas = st.toggle("Show Δm/z", value=False)

        # User selection and plotting
        selected_node_id = selected_feature.split(":")[0]

        # Get all feature IDs in the same cluster as selected one
        selected_cluster = valid_nodes[selected_node_id]
        all_nodes_in_cluster = [
            node_id
            for node_id, cluster in valid_nodes.items()
            if cluster == selected_cluster
        ]

        info = "id" if node_info == "Feature ID" else "mz"
        info_text_col, plot_col = st.columns([1, 4])

        with info_text_col:
            space_for_info = st.empty()
            default_node_colors_dict = {
                "queried_node": "#d2372c",
                "cmmc_match": "#2c921f",
                "fbmn_match": "#c78507",
                "unannotated": "#4b7db4",
            }
            custom_nodes_colors_dict = {}
            with st.expander(":material/palette: Style options"):
                node_size = st.number_input(
                    "Node size (px)", min_value=10, max_value=200, value=60, step=5
                )
                use_custom_node_colors = st.checkbox(
                    "Use custom colors for nodes", key="custom_node_colors"
                )

                for node_type, default_color in default_node_colors_dict.items():
                    node_name = " ".join(node_type.split("_")).upper()
                    custom_nodes_colors_dict[node_type] = st.color_picker(
                        node_name, value=default_color
                    )

                colors_to_use = (
                    custom_nodes_colors_dict if use_custom_node_colors else default_node_colors_dict
                )

            cluster_fig, info_text = plot_cluster_by_node(
                G,
                selected_node_id.split(":")[0],
                all_nodes_in_cluster,
                nodes_info=info,
                node_size=node_size,
                node_colors_dict=colors_to_use,
                show_delta_annotation=show_deltas,
            )
            with space_for_info:
                st.markdown(info_text, unsafe_allow_html=True)

        with plot_col:

            with st.container(border=True):
                st.plotly_chart(cluster_fig.update_layout(dragmode="pan"))

            svg_bytes = cluster_fig.to_image(format="svg")
            st.download_button(
                label=":material/download: Download Plot as SVG",
                data=svg_bytes,
                file_name=f"network_{selected_node_id}.svg",
                mime="image/svg+xml",  # Set the MIME type to SVG
                key='network_plot_download'
            )

        st.markdown("---")

        from microbemass_frame import render_microbemasst_frame
        query_scan_name_mapping = {row["query_scan"]: row["input_name"] for _, row in enriched_result.iterrows()}
        render_microbemasst_frame(input_options_dict=query_scan_name_mapping)

