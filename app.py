import streamlit as st
from streamlit.components.v1 import html

import upset_plot
from box_plot import insert_plot_download_buttons
from network_cluster_plotter import *
from utils import *
from utils import load_uploaded_file_df, validate_task_id_input, fbmn_quant_download_wrapper


def insert_contribute_link(enriched_result, feature_id):
    subset = enriched_result[
        ["LibrarySpectrumID", "query_scan", "input_structure", "input_name", "input_molecule_origin",
         "input_source"]].rename(columns={"LibrarySpectrumID": "input_usi"})
    try:
        params_dict = subset[subset['query_scan'] == int(feature_id.split(":")[0])].to_dict(orient="records")[0]
        params_dict.update({'description': f"Adding information for {feature_id.split(':')[1].strip()}"})
        url = generate_url_hash(params_dict)
        st.markdown(f"- [Contribute depositing more information for {feature_id.split(':')[1]} on CMMC-kb]({url})")
    except IndexError:
        pass


def insert_request_dep_correction_link(enriched_result, feature_id):
    try:
        db_id = enriched_result[enriched_result["query_scan"] == int(feature_id.split(":")[0])]['database_id'].values[0]
        request_correction_subject = urllib.parse.quote(
            f"CMMC-KB Correction request for {feature_id.split(':')[1].strip()} ({db_id})")
        request_correction_body = urllib.parse.quote(
            f"Please provide details about the correction you would like to request for the feature {feature_id.split(':')[1].strip()}\n"
            f"This is the database ID identifier for the deposition you are requesting a correction: {db_id}.\n"
            f"Do not delete it from the email subject.\n\n"
            f"Note that if you just want to include additional information, use the \"Contribute\" link provided in the CMMC dashboard.\n"
            f"This link is only for requesting corrections to the existing data.\n\n"

        )
        st.markdown(
            f"- [Request a correction]"
            f"(mailto:wdnunes@health.ucsd.edu?"
            f"subject={request_correction_subject}"
            f"&body={request_correction_body}"
            f"&cc=hmannochiorusso@health.ucsd.edu)")

    except IndexError:
        pass


def render_details_card(enrich_df, feature_id, columns_to_show):
    """Shows a details card with information about the selected feature."""
    feature_data = enrich_df[enrich_df["query_scan"] == feature_id]
    selected_data = feature_data[columns_to_show]
    try:
        text_info = [
            f"<li><b>{col}</b>: {selected_data.iloc[0][col]}" for col in columns_to_show
        ]
    except IndexError:
        text_info = ["No data available for the selected Feature ID. Probably, the feature ID is not present in the CMMC enrichment results."]
    if not selected_data.empty:
        st.write(f"**Details for Feature ID:** {feature_id}")
        smiles = feature_data.iloc[0]["input_structure"]

        st.image(smiles_to_svg(smiles, (500, 500)))
        st.markdown("<br>".join(text_info), unsafe_allow_html=True)
    else:
        st.warning("No data found for the selected Feature ID.")


def render_color_and_rotation_options(groups, color_prefix, colors_key="overview_plot_custom_check", rotation_key="labels_rot"):
    custom_colors_check = st.checkbox(
        "Use custom colors", key=colors_key
    )
    color_cols = st.columns(3)
    for idx, item in enumerate(groups):
        with color_cols[idx % 3]:
            st.color_picker(
                f"{item}",
                key=f"{color_prefix}_{item}",
                help="Select a color for the group",
                value=st.session_state.get(f"{color_prefix}_{item}", "#1f77b4"),
            )
    rotate_check = st.checkbox('Rotate x-axis labels', value=False, key=f"check_{rotation_key}")
    rotate_labels_angle = st.slider(
        "Rotate x-axis labels",
        min_value=-90,
        max_value=90,
        value=0,
        step=45,
        key=rotation_key,
        label_visibility='collapsed'
    )
    return custom_colors_check, rotate_check, rotate_labels_angle


def render_sidebar():
    global cmmc_task_id, fbmn_task_id, uploaded_quant_file, loaded_metadata_df, include_all_features, run_analysis
    with st.sidebar:
        st.header("📊 Analysis Configuration")
        load_example = st.checkbox("Load Example Data", value=False, key="load_example")
        if not load_example:
            cmmc_task_id = st.text_input(
                "CMMC Enrichment Task ID",
                value=default_cmmc_task_id,
                placeholder="Enter CMMC Enrichment Task ID",
                help="Input your CMMC enrichment task identifier",
                key="cmmc_task_id",
            )
            validate_task_id_input(cmmc_task_id, validation_str="cmmc")

            fbmn_task_id = st.text_input(
                "FBMN Task ID",
                value=default_fbmn_task_id,
                placeholder="Enter FBMN Task ID",
                help="Input your Feature-Based Molecular Network task identifier",
                key="fbmn_task_id",
            )
            validate_task_id_input(fbmn_task_id, 'feature_based')

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
            if uploaded_metadata_file is not None:
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

                    # Show preview
                    with st.expander("Preview Data", icon=":material/visibility:"):
                        st.dataframe(loaded_metadata_df.head(), use_container_width=True)

                except Exception as e:
                    st.error(f"Error reading file: {str(e)}", icon=":material/error:")
            else:
                st.info("📤 Please upload a metadata table")
        else:
            # loads Quinn's 2020 example data https://doi.org/10.1038/s41586-020-2047-9
            cmmc_task_id = "7f53b63490c945e980dfa10273a296cd"
            fbmn_task_id = "58e0e2959ec748049cb2c5f8bb8b87dc"
            st.session_state['fbmn_task_id'] = fbmn_task_id
            st.session_state['cmmc_task_id'] = cmmc_task_id
            uploaded_metadata_file = open('data/metadata_quinn2020.tsv', 'rb')

            st.write("Using example data from Quinn et al. 2020: https://doi.org/10.1038/s41586-020-2047-9")
            st.write("CMMC Task ID: ", cmmc_task_id)
            st.write("FBMN Task ID: ", fbmn_task_id)

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


# TODO: Bump version
app_version = "2025-06-30"

menu_items={"about": f"**CMMC Dashboard version: {app_version}**"}
st.set_page_config(
    page_title="CMMC Analysis Dashboard", page_icon="favicon.png", layout="wide", menu_items=menu_items
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
        merged_df = box_plot.prepare_lcms_data(
            df_quant, loaded_metadata_df, enriched_result, include_all_features
        )


        progress_bar.progress(85)

        # Step 5: Fetch network data
        status_text.text("Fetching molecular network data...")
        graphml_file_name = fetch_cmmc_graphml(
            cmmc_task_id
        )

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

    render_welcome_message()

# Main content area
if st.session_state.get("run_analysis"):
    st.title("🦠 CMMC Analysis Dashboard")
    tabs = st.tabs(["Data Explorer", 'Advanced Visualizations'])
    with tabs[0]:

        st.subheader(
            "Data Overview Box Plots",
            help="Select the column that contains the groups you want to compare to see a boxplot for each detected feature.",
        )
        data_overview_df = st.session_state.get("merged_df")

        col1, col2 = st.columns(2)
        with col1:
            column_select = st.selectbox(
                "Select column", sorted([i for i in data_overview_df.columns])
            )
            if column_select:
                group_by = st.multiselect(
                    "Select groups to compare",
                    [i for i in data_overview_df[column_select].unique()],
                    key="a",
                )
        with col2:
            st.write(
                '<div style="height: 25px;"></div>', unsafe_allow_html=True
            )  # this is just to align the button with the textinput field
            # Filter data_overview_df based on the selected column and value
            with st.expander(
                    "Filter options - Source or Origin", icon=":material/filter_alt:"
            ):
                input1, input2 = st.columns(2)
                with input1:
                    first = st.selectbox(
                        "Column", ["input_molecule_origin", "input_source"]
                    )
                with input2:
                    origin_list = [
                        "Ambiguous",
                        "De novo biosynthesis by microbes",
                        "Diet",
                        "Drug",
                        "Exposure",
                        "Exposure/diet",
                        "Host",
                        "Host metabolism of microbial metabolites",
                        "Insecticides/pesticides",
                        "Microbial metabolism of drugs",
                        "Microbial metabolism of food molecules",
                        "Microbial metabolism of host-derived molecules",
                        "Microbial metabolism of microbial-derived molecules",
                        "Microbial metabolism of other human-made molecules",
                        "Unknown/Undefined",
                    ]
                    source_list = [
                        "Microbial",
                        "Host",
                        "Diet",
                        "Unknown",
                        "Ambiguous",
                        "Drug",
                        "Exposure",
                        "Pesticides/insecticides",
                        "Other human-made molecules",
                    ]
                    second = st.multiselect(
                        "Value",
                        (
                            origin_list
                            if first == "input_molecule_origin"
                            else source_list
                        ),
                    )
                st.write(
                    "**Tip**: use the [UpSet Plot](#up-set-plot) to see the possible groupings"
                )
                from utils import prepare_dataframe, find_exact_matches

                filter_results = render_filter_options(
                    data_overview_df, first, second, key="overview"
                )
            st.write('<div style="height: 20px;"></div>', unsafe_allow_html=True)
            with st.expander("Style Options", icon=":material/palette:"):

                groups = [i for i in group_by]
                custom_colors, rotate_check, rotate_labels_angle = render_color_and_rotation_options(groups, "color",
                                                                                rotation_key="overview_labels_rot")
        data_overview_df = filter_results.data
        filter_string = filter_results.filters

        feat_id_dict = dict(
            zip(data_overview_df["featureID"], data_overview_df["input_name"])
        )
        feat_id_dict = {str(k): str(v) for k, v in feat_id_dict.items()}
        feat_id_dict = dict(sorted(feat_id_dict.items(), key=lambda item: item[1]))

        fid_items = [f"{k}: {v}" for k, v in feat_id_dict.items()]
        col_fid_1, col_download = st.columns([3, 1])
        with col_fid_1:
            feature_id = st.selectbox(
                f"Select Feature ID :blue-badge[{len(fid_items)} item(s)]",
                [None] + fid_items,
                key="b",
            )

        with col_download:
            group_colors = {
                group: st.session_state.get(f"color_{group}", "#1f77b4")
                for group in group_by
            }

            if len(data_overview_df) > 0 and len(feat_id_dict) > 0:
                st.write(
                    '<div style="height: 28px;"></div>', unsafe_allow_html=True
                )  # this is just to align the button with the textinput field
                add_pdf_download_overview(
                    data_overview_df,
                    feat_id_dict,
                    group_by,
                    column_select,
                    filter_string,
                    color_mapping=group_colors if custom_colors else None
                )
        if len(data_overview_df) > 0:
            if feature_id:
                plot_col, details_col = st.columns([3, 1])
                with plot_col:
                    overview_plot, overview_plot_df = box_plot.plot_boxplots_by_group(
                        data_overview_df,
                        groups1=group_by,  # this will be on x axis
                        column1=column_select,
                        feature_id=int(feature_id.split(":")[0]),
                        informations=filter_string,
                        color_mapping=group_colors if custom_colors else None,
                    )

                    if rotate_check:
                        overview_plot.update_xaxes(tickangle=rotate_labels_angle)
                    st.plotly_chart(
                        overview_plot,
                        use_container_width=True,
                        key="graph1",
                    )
                    svg_bytes = overview_plot.to_image(format="svg")
                    #add two buttons for downloading plot figure and data
                    insert_plot_download_buttons(overview_plot_df, feature_id, svg_bytes, key_prefix="overview")

                with details_col:
                    # Show details card for the selected feature ID
                    enriched_result = st.session_state.get("enriched_result")
                    with st.expander("Details", icon=":material/info:"):
                        columns_to_show = st.multiselect(
                            "Select columns to show in details card",
                            enriched_result.columns.tolist(),
                            default=[
                                "input_name",
                                "input_molecule_origin",
                                "input_source",
                            ],
                        )
                    render_details_card(
                        enriched_result, int(feature_id.split(":")[0]), columns_to_show
                    )
                insert_contribute_link(enriched_result, feature_id)

                # renders a link to request a correction
                insert_request_dep_correction_link(enriched_result, feature_id)

            else:
                st.warning("Select a feature ID to plot", icon="🆔")
        else:
            st.warning(
                "The selected filter did not return any results. Try again with another combination"
            )

        if st.session_state.get("run_analysis"):
            st.markdown("---")
            st.subheader(
                "📊 Stratified Box Plots",
                help="**Group 1:** Stratify the data for the selected attribute. **Group 2:** Select the groups to visualize",
            )

            metadata = st.session_state.get("metadata_df")
            merged_data = st.session_state.get("merged_df")

            col_attr1, col_attr2 = st.columns(2)
            with col_attr1:
                selected_attribute1 = st.selectbox(
                    "Metadata prefilter (group 1)",
                    [None] + sorted([i for i in metadata.columns]),
                    help="**Group 1:** Stratify the data for the selected attribute",
                )
            with col_attr2:
                selected_attribute2 = st.selectbox(
                    "Metadata group 2 (x axis)",
                    sorted([i for i in metadata.columns]),
                    help="**Group 2:** Select the groups to visualize",
                )

            col1, col2 = st.columns([1, 1])
            with col1:
                if selected_attribute1:
                    groups1 = st.selectbox(
                        "Group 1 (single selection)",
                        [i for i in metadata[selected_attribute1].unique()],
                    )
                else:
                    groups1 = None

            with col2:
                if selected_attribute2:
                    groups2 = st.multiselect(
                        "Group 2",
                        [i for i in metadata[selected_attribute2].unique()],
                        placeholder="Choose one or more",
                    )
            filter_col, style_col = st.columns([1, 1])
            with filter_col:
                with st.expander(
                        "Filter options - Source or Origin", icon=":material/filter_alt:"
                ):
                    input1, input2 = st.columns(2)
                    with input1:
                        first = st.selectbox(
                            "Column",
                            ["input_molecule_origin", "input_source"],
                            key="firstB",
                        )
                    with input2:
                        origin_list = [
                            "Ambiguous",
                            "De novo biosynthesis by microbes",
                            "Diet",
                            "Drug",
                            "Exposure",
                            "Exposure/diet",
                            "Host",
                            "Host metabolism of microbial metabolites",
                            "Insecticides/pesticides",
                            "Microbial metabolism of drugs",
                            "Microbial metabolism of food molecules",
                            "Microbial metabolism of host-derived molecules",
                            "Microbial metabolism of microbial-derived molecules",
                            "Microbial metabolism of other human-made molecules",
                            "Unknown/Undefined",
                        ]
                        source_list = [
                            "Microbial",
                            "Host",
                            "Diet",
                            "Unknown",
                            "Ambiguous",
                            "Drug",
                            "Exposure",
                            "Pesticides/insecticides",
                            "Other human-made molecules",
                        ]
                        second = st.multiselect(
                            "Value",
                            (
                                origin_list
                                if first == "input_molecule_origin"
                                else source_list
                            ),
                            key="secondB",
                        )

                    filtered_results_boxplot = render_filter_options(
                        merged_data, first, second, key="Boxplots"
                    )
                    merged_data = filtered_results_boxplot.data
                    boxp_filter_string = filtered_results_boxplot.filters

            with style_col:
                with st.expander("Style Options", icon=":material/palette:"):
                    groups = [i for i in groups2]
                    custom_colors, rotate_check, rotate_labels_angle = render_color_and_rotation_options(groups, "box_color",
                                                                                    rotation_key="boxp_labels_rot", colors_key="box_plot_custom_check")

            # from merged data create an input widget to select featureID (with input_name) from merged_data.columns
            feat_id_dict = dict(zip(merged_data["featureID"], merged_data["input_name"]))
            feat_id_dict = {str(k): str(v) for k, v in feat_id_dict.items()}
            feat_id_dict = dict(sorted(feat_id_dict.items(), key=lambda item: item[1]))

            prefilter = selected_attribute1 if selected_attribute1 != "None" else None
            if prefilter:
                boxp_filter_string += f" | prefilter: {prefilter} = {groups1}"

            col_fid, col_download = st.columns([3, 1])

            with col_fid:
                fid_items_2 = [f"{k}: {v}" for k, v in feat_id_dict.items()]
                feature_id = st.selectbox(
                    f"Select Feature ID :blue-badge[{len(fid_items_2)} item(s)]",
                    [None] + fid_items_2,
                )

            # Add PDF download button
            with col_download:
                group_colors = {
                    group: st.session_state.get(f"box_color_{group}", "#1f77b4")
                    for group in groups2
                }
                st.write(
                    '<div style="height: 28px;"></div>', unsafe_allow_html=True
                )  # this is just to align the button with the textinput field
                add_pdf_download_boxplots(
                    merged_data,
                    feat_id_dict,
                    groups1,
                    groups2,
                    selected_attribute2,
                    prefilter,
                    str(boxp_filter_string),
                    color_mapping=group_colors if custom_colors else None,
                )

            if feature_id:
                plot_col, details_col = st.columns([3, 1])
                with plot_col:
                    boxplot_plot, box_plot_data = box_plot.plot_boxplots_by_group(
                        merged_data,
                        groups2,
                        [groups1],
                        int(feature_id.split(":")[0]),
                        selected_attribute2,
                        prefilter,
                        informations=boxp_filter_string,
                        color_mapping=group_colors if custom_colors else None,
                    )
                    if rotate_check:
                        boxplot_plot.update_xaxes(tickangle=rotate_labels_angle)
                    st.plotly_chart(
                        boxplot_plot,
                        use_container_width=True,
                        key="graph2",
                    )
                    svg_bytes = boxplot_plot.to_image(format="svg")
                    insert_plot_download_buttons(box_plot_data, feature_id, svg_bytes, key_prefix="stratified")
                with details_col:
                    # Show details card for the selected feature ID
                    enriched_result = st.session_state.get("enriched_result")
                    with st.expander("Details", icon=":material/info:"):
                        columns_to_show = st.multiselect(
                            "Select columns to show in details card",
                            enriched_result.columns.tolist(),
                            default=["input_name", "input_molecule_origin", "input_source"],
                            key="details_box_plot_columns",
                        )
                    render_details_card(
                        enriched_result, int(feature_id.split(":")[0]), columns_to_show
                    )

                insert_contribute_link(enriched_result, feature_id)

                #renders a link to request a correction
                insert_request_dep_correction_link(enriched_result, feature_id)

            else:
                st.warning("Select all required fields to see the boxplot")

        if st.session_state.get("run_analysis"):
            st.markdown("---")
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
                _, plot_col, _ = st.columns([1, 1, 1])
                upset_fig = upset_fig_source
            else:
                _, plot_col, _ = st.columns([1, 4, 1])
                upset_fig = upset_fig_origin

            with plot_col:
                st.image(upset_fig, use_container_width=False)
                st.download_button(
                    label=":material/download: Download as SVG",
                    data=upset_fig,
                    file_name="upset_plot.svg",
                    mime="image/svg+xml",
                    key='upset_plot_download'
                )

    # INTERFACE ELEMENTS
    st.markdown("---")
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

    with tabs[1]:
        st.subheader("🕸️ Molecular Network Visualization")

        col_select, col_radio, col_deltas  = st.columns([1, 1, 1])
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
            show_deltas = st.checkbox("Show Δm/z", value=False)


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

        from microbemass_frame import render_microbemasst_frame
        st.markdown("---")
        render_microbemasst_frame()

