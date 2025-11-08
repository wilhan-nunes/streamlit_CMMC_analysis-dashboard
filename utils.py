import base64
import gzip
import json
import os
import tempfile
import urllib.parse
from dataclasses import dataclass
from io import BytesIO
from typing import Literal

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import plotly.io as pio
import requests
import streamlit as st
from gnpsdata import taskinfo, taskresult, workflow_fbmn
from matplotlib.backends.backend_pdf import PdfPages
from rdkit import Chem
from rdkit.Chem import Draw

import box_plot


def get_git_short_rev():
    try:
        with open('.git/logs/HEAD', 'r') as f:
            last_line = f.readlines()[-1]
            hash_val = last_line.split()[1]
        return hash_val[:7]
    except Exception:
        return ".git/ not found"


def generate_url_hash(params_dict):
    url = "https://gnps2.org/workflowinput?workflowname=cmmc_deposition_workflow"
    url += "#{}".format(urllib.parse.quote(json.dumps(params_dict)))
    return url


def smiles_to_structure_image(smiles_string, img_size=(300, 300), save_path=None):
    """
    Convert a SMILES string to a molecular structure image.

    Args:
        smiles_string (str): The SMILES representation of the molecule
        img_size (tuple): Size of the output image (width, height)
        save_path (str, optional): Path to save the image file

    Returns:
        PIL.Image: The molecular structure image
    """
    try:
        # Parse the SMILES string into a molecule object
        mol = Chem.MolFromSmiles(smiles_string)

        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles_string}")

        # Generate 2D coordinates for the molecule
        Chem.rdDepictor.Compute2DCoords(mol)

        # Create the molecular structure image
        img = Draw.MolToImage(mol, size=img_size)

        # Save the image if a path is provided
        if save_path:
            img.save(save_path)
            print(f"Image saved to: {save_path}")

        return img

    except Exception as e:
        print(f"Error converting SMILES to image: {e}")
        return None


def smiles_to_inchikey(smiles_string):
    """
    Convert a SMILES string to an InChIKey.

    Args:
        smiles_string (str): The SMILES representation of the molecule

    Returns:
        str: The InChIKey of the molecule
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)

        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles_string}")

        inchikey = Chem.MolToInchiKey(mol)
        return inchikey

    except Exception as e:
        print(f"Error converting SMILES to InChIKey: {e}")
        return None


def smiles_to_svg(smiles_string, svg_size=(300, 300)):
    """
    Convert a SMILES string to an SVG molecular structure.

    Args:
        smiles_string (str): The SMILES representation of the molecule
        svg_size (tuple): Size of the SVG (width, height)

    Returns:
        str: SVG string representation of the molecule
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)

        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles_string}")

        Chem.rdDepictor.Compute2DCoords(mol)

        # Generate SVG
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(svg_size[0], svg_size[1])
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        return svg

    except Exception as e:
        print(f"Error converting SMILES to SVG: {e}")
        return None



@st.cache_data(show_spinner=False)
def fetch_enriched_results(task_id: str) -> pd.DataFrame:
    """
    Download enriched results as a tsv file from GNPS2 Enrichment workflow.

    :param task_id: GNPS2 Enrichment workflow task ID
    :return: pd.DataFrame containing enriched results
    """
    url = taskresult.determine_gnps2_resultfile_url(task_id, 'nf_output/cmmc_results/cmmc_enriched_results.tsv')
    df = pd.read_csv(url, sep="\t", low_memory=False)
    return df


@st.cache_data(show_spinner=False)
def fetch_phylogeny_results(task_id: str) -> pd.DataFrame:
    """
    Fetch phylogeny results from GNPS2.

    :param task_id: GNPS2 Enrichment workflow task ID
    :return: pd.DataFrame containing phylogeny results
    """
    url = taskresult.determine_gnps2_resultfile_url(task_id, "nf_output/cmmc_results/cmmc_taxonomy.tsv")
    df = pd.read_csv(url, sep="\t", low_memory=False)
    return df


@st.cache_data(show_spinner=False)
def fetch_cmmc_graphml(task_id: str):
    url = taskresult.determine_gnps2_resultfile_url(task_id, "nf_output/gnps_network/network.graphml")
    response = requests.get(url)
    if response.status_code == 200:
        return BytesIO(response.content)
    else:
        raise Exception(f"Failed to fetch graphml file: {response.status_code}")


# this is used for the FBMN files
@st.cache_data(show_spinner=False)
def fbmn_quant_download_wrapper(task_id):
    return workflow_fbmn.get_quantification_dataframe(task_id, gnps2=True)


@st.cache_data(show_spinner=False)
def prepare_dataframe(enrich_df, by: Literal["source", "origin"]):
    """Prepare dataframe for filtering metabolites sources and origins.
    This could be then used to filter the original dataframe"""
    # Select the column1 to process
    if by == "source":
        column = "input_source"
    elif by == "origin":
        column = "input_molecule_origin"
    else:
        raise ValueError("Parameter 'by' must be either 'source' or 'origin'.")

    # Clean and standardize the selected column1, ensuring unique values per row
    enrich_df["input_clean"] = (
        enrich_df[column]
        .fillna("")
        .str.replace(r"\s+and\s+", ";", regex=True)
        .str.split(";")
        .apply(lambda items: list({item.strip() for item in items if item}))
    )
    all_categories = sorted(
        {item.strip() for items in enrich_df["input_clean"] for item in items if item}
    )

    df_indicators = pd.DataFrame()
    for category in all_categories:
        df_indicators[category] = enrich_df["input_clean"].apply(
            lambda x: category in x
        )

    df_indicators = df_indicators[df_indicators.any(axis=1)]
    return df_indicators


@st.cache_data(show_spinner=False)
def find_exact_matches(df, target_cols):
    """
    Function tah receives the enrichment results dataframe and filter the rows according to the target_cols
    :param df: enrichment results dataframe
    :param target_cols: columns from the prepared dataframe containing True/False rows for each category (column)
    :return: pd.DataFrame with booleans indicating to which categories each row pertains.
    """

    try:
        target_mask = (df[target_cols] == True).all(axis=1)
        other_cols = [col for col in df.columns if col not in target_cols]
        other_mask = (df[other_cols] == False).all(axis=1)
        exact_match = target_mask & other_mask
        return df[exact_match]
    except KeyError:
        # return empty dataframe in case of no matches
        return pd.DataFrame()


def load_uploaded_file_df(uploaded_file):
    if uploaded_file.name.endswith(".csv"):
        loaded_file_df = pd.read_csv(uploaded_file)
    elif uploaded_file.name.endswith(".tsv") or uploaded_file.name.endswith(".txt"):
        loaded_file_df = pd.read_csv(uploaded_file, sep="\t")
    else:  # Excel files
        loaded_file_df = pd.read_excel(uploaded_file)
    return loaded_file_df


def validate_task_id_input(task_id: str, validation_str: str):
    if task_id and len(task_id) == 32:
        try:
            task_data = taskinfo.get_task_information(task_id)
            workflow_name = task_data['workflowname']
            if "enrichment" in workflow_name:
                st.session_state['enrichment_date'] = task_data['create_time'].split(' ')[0]
            if validation_str not in workflow_name:
                st.error(f"The provided task ID is for {workflow_name}", icon=":material/error:")
        except Exception:
            st.error(f"Failed to fetch task information. Is this a valid task id?", icon=":material/dangerous:")
    elif task_id:
        st.warning("Task ID must be exactly 32 characters long.")


@st.cache_data(show_spinner=False)
def prepare_lcms_data(
        df_quant: pd.DataFrame, df_metadata: pd.DataFrame, cmmc_results: pd.DataFrame, include_all_scans: bool = False
):
    """
    Optimized version of prepare_lcms_data with improved performance for large datasets.
    """
    # Work with a copy to avoid modifying original data
    df_quant = df_quant.copy()

    # Early filtering to reduce dataset size
    if not include_all_scans:
        microbial_scans = set(cmmc_results["query_scan"].tolist())  # Use set for faster lookup
        df_quant = df_quant[df_quant["row ID"].isin(microbial_scans)]

    # Optimize column operations - do them all at once
    df_quant.columns = df_quant.columns.str.replace(" Peak area", "", regex=False)

    # More efficient transpose and filtering
    # First, identify the data rows we need (those containing mzML/mzXML/.raw)
    data_mask = df_quant.columns.str.contains("mzML|mzXML|\\.raw", case=False, na=False)
    if not data_mask.any():
        # If no columns match, check the index after transpose
        t_df_quant = df_quant.set_index("row ID").T
        data_mask = t_df_quant.index.str.contains("mzML|mzXML|\\.raw", case=False, na=False)
        t_df_quant = t_df_quant[data_mask].reset_index()
        t_df_quant.columns.name = None
        t_df_quant = t_df_quant.rename(columns={"index": "filename"})
    else:
        # Standard transpose operation
        t_df_quant = df_quant.set_index("row ID").T
        data_mask = t_df_quant.index.str.contains("mzML|mzXML|\\.raw", case=False, na=False)
        t_df_quant = t_df_quant[data_mask].reset_index()
        t_df_quant.columns.name = None
        t_df_quant = t_df_quant.rename(columns={"index": "filename"})

    # Optimize the melting operation
    # Use only numeric columns for melting (exclude filename)
    numeric_cols = [col for col in t_df_quant.columns if col != "filename"]

    # More efficient melt with explicit column specification
    df_quant_long = pd.melt(
        t_df_quant,
        id_vars=["filename"],
        value_vars=numeric_cols,
        var_name="featureID",
        value_name="Peak Area"
    )

    # Convert featureID to numeric if it isn't already (for better merge performance)
    df_quant_long["featureID"] = pd.to_numeric(df_quant_long["featureID"], errors='coerce')

    # Optimize filename cleaning - do both at once with more efficient regex
    filename_pattern = r"\.mzML|\.mzXML|\.raw"
    df_metadata = df_metadata.copy()
    df_metadata['filename'] = df_metadata['filename'].str.replace(filename_pattern, "", regex=True)
    df_quant_long['filename'] = df_quant_long['filename'].str.replace(filename_pattern, "", regex=True)

    # Perform merges with optimized settings
    # First merge: metadata with quantification data
    df_quant_merged = pd.merge(
        df_quant_long,
        df_metadata,
        on="filename",
        how="left"  # Use inner join to reduce size if possible
    )

    # Optimize the CMMC results merge
    cmmc_results_clean = cmmc_results.rename(columns={"query_scan": "featureID"})

    # Use left merge but optimize by ensuring data types match
    if not include_all_scans:
        # When not including all scans, we can use inner join for better performance
        df_quant_merged = pd.merge(
            df_quant_merged,
            cmmc_results_clean,
            on="featureID",
            how="inner"
        )
    else:
        # When including all scans, use left join but optimize column selection
        # Only keep essential columns from cmmc_results to reduce memory usage
        essential_cmmc_cols = ["featureID", "input_name", "input_molecule_origin", "input_source"]
        available_cols = [col for col in essential_cmmc_cols if col in cmmc_results_clean.columns]

        df_quant_merged = pd.merge(
            df_quant_merged,
            cmmc_results_clean[available_cols],
            on="featureID",
            how="left"
        )

    # Optimize data types to reduce memory usage
    # Convert categorical columns to category dtype
    # categorical_columns = ['input_molecule_origin', 'input_source', 'input_name']
    # for col in categorical_columns:
    #     if col in df_quant_merged.columns:
    #         df_quant_merged[col] = df_quant_merged[col].astype('category')

    # Convert abundance to float32 if it's float64 (reduces memory by ~50%)
    if df_quant_merged['Peak Area'].dtype == 'float64':
        df_quant_merged['Peak Area'] = df_quant_merged['Peak Area'].astype('float32')

    return df_quant_merged


def insert_contribute_link(enriched_result, feature_id):
    subset = enriched_result[
        ["LibrarySpectrumID", "query_scan", "input_structure", "input_name", "input_molecule_origin",
         "input_source"]].rename(columns={"LibrarySpectrumID": "input_usi"})
    try:
        params_dict = subset[subset['query_scan'] == int(feature_id.split(":")[0])].to_dict(orient="records")[0]
        params_dict.update({'description': f"Adding information for {feature_id.split(':')[1].strip()}"})
        url = generate_url_hash(params_dict)
        st.link_button(
            label=f"Contribute info for {feature_id.split(':')[1]}",
            url=url,
            icon=":material/database_upload:"
        )
    except IndexError:
        pass


def insert_request_dep_correction_link(enriched_result, feature_id):
    try:
        db_id = enriched_result[enriched_result["query_scan"] == int(feature_id.split(":")[0])]['database_id'].values[0]
        request_correction_subject = urllib.parse.quote(
            f"CMMC-KB Correction request for {feature_id.split(':')[1].strip()} ({db_id})")
        request_correction_body = urllib.parse.quote(
            f"Please, read the instructions below. You can then clean the body of the email to include your request details.\n\n"
            f"Please provide details about the correction you would like to request for the feature {feature_id.split(':')[1].strip()}\n"
            f"This is the database ID identifier for the deposition you are requesting a correction: {db_id}.\n"
            f"Do not delete it from the email subject.\n\n"
            f"Note that if you just want to include additional information, use the \"Contribute\" link provided in the CMMC dashboard.\n"
            f"This link is only for requesting corrections to the existing data.\n\n"
            f"Thank you for your contribution to the CMMC knowledge base!\n\n"

        )
        st.link_button(
            label="Request a correction",
            url=f"mailto:wdnunes@health.ucsd.edu?subject={request_correction_subject}&body={request_correction_body}&cc=hmannochiorusso@health.ucsd.edu",
            icon=":material/ink_eraser:")

    except IndexError:
        pass


def render_details_card(enrich_df, feature_id, columns_to_show, cmmc_task_id):
    """Shows a details card with information about the selected feature."""
    feature_data = enrich_df[enrich_df["query_scan"] == feature_id]
    selected_data = feature_data[columns_to_show]
    try:
        text_info = []
        for col in columns_to_show:
            if "clean" in col:
                # Join tuple/list items with semicolon
                items_str = "; ".join(list(selected_data.iloc[0][col]))
                col_display = col.replace('_clean', '')
                text_info.append(f"- **{col_display}**: {items_str}")
            else:
                text_info.append(f"- **{col}**: {selected_data.iloc[0][col]}")

    except IndexError:
        text_info = ["No data available for the selected Feature ID. Probably, the feature ID is not present in the CMMC enrichment results."]
    if not selected_data.empty:
        lib_usi = feature_data.iloc[0]["LibrarySpectrumID"]
        task_usi = f"mzspec:GNPS2:TASK-{cmmc_task_id}-nf_output/gnps_network/specs_ms.mgf:scan:{feature_id}"
        base_url = "https://metabolomics-usi.gnps2.org/dashinterface/"
        mirror_plot_link = base_url + f"?usi1={lib_usi}&usi2={task_usi}"

        st.write(f"**Details for Feature ID:** {feature_id} - [Mirror plot view]({mirror_plot_link})")
        smiles = feature_data.iloc[0]["input_structure"]
        inchikey = smiles_to_inchikey(smiles)
        enrichment_date_ = st.session_state.get('enrichment_date', 'N/A')
        smiles_svg = smiles_to_svg(smiles, (500, 500))
        latest_data_info = f"""Data as of **{enrichment_date_}** - [View latest data](https://cmmc-kb.gnps2.org/structurepage/?inchikey={inchikey})"""

        if smiles_svg:
            with st.container(border=True, key="structure_image"):
                st.image(smiles_svg)
            st.info(latest_data_info)
        else:
            st.info("No valid SMILES string available to generate structure image.")
            st.info(latest_data_info)
    
        st.markdown("\n".join(text_info))
    else:
        st.warning("No data found for the selected Feature ID.")



def get_session_state_params():
    """
    Retrieve all session state parameters excluding data, dataframes, and network objects
    """
    # Define keys to exclude (data-related items)
    exclude_keys = {
        'enriched_result', 'df_quant', 'merged_df', 'G', 'metadata_df',
        # Add any other data-related keys you want to exclude
        'uploaded_file', 'loaded_data', 'network_data', 'network_plot_download'
    }

    # Get all session state items except excluded ones
    params = {
        key: value for key, value in st.session_state.items()
        if key not in exclude_keys and not isinstance(value, (pd.DataFrame, nx.Graph))
    }

    return params

def generate_boxplot_script(
        feature_id,
        grouping_column,
        selected_groups,
        intensity_col,
        stratify_column,
        selected_strata,
        selected_test,
        alpha_level,
        use_custom_colors,
        custom_colors,
        use_log_scale,
        rotate_angle,
        show_grid=True
):
    """
    Generate a complete Python script for recreating the app boxplots.

    Args:
        feature_id: ID of the feature being plotted
        grouping_column: Column name for grouping
        selected_groups: List of selected groups
        intensity_col: Column name for intensity values
        stratify_column: Column name for stratification (can be None)
        selected_strata: List of selected strata (can be None)
        selected_test: Statistical test to perform
        alpha_level: Alpha level for significance testing
        use_custom_colors: Boolean for custom color usage
        custom_colors: Dictionary of custom colors
        use_log_scale: Boolean for log scale usage
        rotate_angle: Angle for x-axis label rotation
        show_grid: Boolean for showing grid lines

    Returns:
        str: Complete Python script as a string
    """
    file_suffix = 'paired' if stratify_column else 'simple'

    script = f"""import pandas as pd
import plotly.express as px
from scipy.stats import mannwhitneyu, kruskal, ttest_ind, f_oneway

# Load CSV data -> use either the csv file exported for a single feature or the full CMMC enrichment data file
input_file = 'filtered_data_{grouping_column}_{file_suffix}.csv'
df = pd.read_csv(input_file, sep=',')

# Filter for feature and groups
feature_id = {repr(feature_id)}
grouping_column = {repr(grouping_column)}
selected_groups = {repr(selected_groups)}
intensity_col = {repr(intensity_col)}
stratify_column = {repr(stratify_column)}
# Set selected_strata to None or [] if not stratifying
selected_strata = {repr(selected_strata)}
selected_test = {repr(selected_test)}
alpha_level = {repr(alpha_level)}

# Apply additional grouping and stratification filters
filter_conditions = [
    (df['featureID'] == feature_id),
    (df[grouping_column].isin(selected_groups))
]

plot_data = df[
    filter_conditions[0] & filter_conditions[1] &
    (filter_conditions[2] if len(filter_conditions) > 2 else True)
    ].copy()

if stratify_column and selected_strata:
    plot_data = plot_data[plot_data[stratify_column].isin(selected_strata)]

# Style options
use_custom_colors = {repr(use_custom_colors)}
custom_colors = {repr(custom_colors)}
use_log_scale = {repr(use_log_scale)}
rotate_angle = {repr(rotate_angle)}
show_grid = {repr(show_grid)}

# Apply log scale transformation if needed
if use_log_scale:
    plot_data[intensity_col] = plot_data[intensity_col].apply(
        lambda x: x if x > 0 else 1e-9  # Avoid log(0)
    )

# Category orders for consistent group/strata ordering
category_orders = {{grouping_column: selected_groups}}
if stratify_column and selected_strata:
    category_orders[stratify_column] = selected_strata

# Statistical test function
def run_statistical_test(data, grouping_column, intensity_col, selected_groups, test_type, alpha, stratify_column=None, selected_strata=None):
    results = dict()
    def get_significance(p):
        return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    if stratify_column and selected_strata:
        results = dict()
        for stratum in selected_strata:
            stratum_data = data[data[stratify_column] == stratum]
            group_data = [stratum_data[stratum_data[grouping_column] == g][intensity_col].dropna() for g in selected_groups]
            group_data = [g for g in group_data if len(g) > 0]
            if len(group_data) < 2:
                results[stratum] = {{'error': f'Need at least 2 groups with data for stratum: {{stratum}}'}}
                continue
            if test_type == 'Mann-Whitney U (2 groups)':
                stat, p = mannwhitneyu(group_data[0], group_data[1], alternative='two-sided')
                test_name = 'Mann-Whitney U'
            elif test_type == 'Kruskal-Wallis (>= 3 groups)':
                stat, p = kruskal(*group_data)
                test_name = 'Kruskal-Wallis'
            elif test_type == 'T-test (2 groups)':
                stat, p = ttest_ind(group_data[0], group_data[1])
                test_name = 'T-test'
            elif test_type == 'ANOVA (>= 3 groups)':
                stat, p = f_oneway(*group_data)
                test_name = 'ANOVA'
            else:
                stat, p, test_name = None, None, 'Unknown'
            results[stratum] = {{
                'test': test_name,
                'statistic': stat,
                'p_value': p,
                'significance': get_significance(p)
            }}
    else:
        group_data = [data[data[grouping_column] == g][intensity_col].dropna() for g in selected_groups]
        group_data = [g for g in group_data if len(g) > 0]
        if len(group_data) < 2:
            return {{'error': 'Need at least 2 groups with data for statistical testing'}}
        if test_type == 'Mann-Whitney U (2 groups)':
            stat, p = mannwhitneyu(group_data[0], group_data[1], alternative='two-sided')
            test_name = 'Mann-Whitney U'
        elif test_type == 'Kruskal-Wallis (>= 3 groups)':
            stat, p = kruskal(*group_data)
            test_name = 'Kruskal-Wallis'
        elif test_type == 'T-test (2 groups)':
            stat, p = ttest_ind(group_data[0], group_data[1])
            test_name = 'T-test'
        elif test_type == 'ANOVA (>= 3 groups)':
            stat, p = f_oneway(*group_data)
            test_name = 'ANOVA'
        else:
            stat, p, test_name = None, None, 'Unknown'
        results = {{
            'test': test_name,
            'statistic': stat,
            'p_value': p,
            'significance': get_significance(p)
        }}
    return results

# Create boxplot with consistent styles
fig = px.box(
    plot_data,
    x=grouping_column,
    y=intensity_col,
    color=grouping_column,
    facet_col=stratify_column if stratify_column else None,
    points="all",
    color_discrete_map=custom_colors if use_custom_colors else None,
    color_discrete_sequence=None if use_custom_colors else px.colors.qualitative.Set1,
    category_orders=category_orders
)
fig.update_layout(
    title=f"Boxplot for Feature {{feature_id}}: {{plot_data[plot_data['featureID'] == feature_id]['input_name'].iloc[0]}}",
    xaxis_title=grouping_column,
    yaxis_title=intensity_col,
    template="plotly_white",
    showlegend=False,
    height=600,
    yaxis=dict(showgrid=show_grid)
)
fig.update_xaxes(tickangle=rotate_angle)
if use_log_scale:
    fig.update_yaxes(type="log", exponentformat="power", showexponent="all")
fig.update_traces(
    boxpoints='all',
    marker=dict(opacity=0.6),
    pointpos=0,
    jitter=0.3,
)

# Statistical annotation
results = run_statistical_test(plot_data, grouping_column, intensity_col, selected_groups, selected_test, alpha_level, stratify_column, selected_strata)
if stratify_column and selected_strata:
    for i, stratum in enumerate(selected_strata, start=1):
        if stratum in results and 'error' not in results[stratum]:
            stat = results[stratum]['statistic']
            p = results[stratum]['p_value']
            sig = results[stratum]['significance']
            y_max = plot_data[plot_data[stratify_column] == stratum][intensity_col].max()
            y_min = plot_data[plot_data[stratify_column] == stratum][intensity_col].min()
            y_range = y_max - y_min
            fig.add_annotation(
                x=0.5, y=y_max + y_range * 0.1,
                xref="paper", yref="y",
                text=f"{{results[stratum]['test']}}: p = {{p:.4f}} ({{sig}})",
                showarrow=False,
                font=dict(size=12),
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="black",
                borderwidth=1,
                row=1, col=i
            )
            fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
else:
    if 'error' not in results:
        stat = results['statistic']
        p = results['p_value']
        sig = results['significance']
        y_max = plot_data[intensity_col].max()
        y_min = plot_data[intensity_col].min()
        y_range = y_max - y_min
        fig.add_annotation(
            x=0.5, y=y_max + y_range * 0.1,
            xref="paper", yref="y",
            text=f"{{results['test']}}: p = {{p:.4f}} ({{sig}})",
            showarrow=False,
            font=dict(size=12),
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="black",
            borderwidth=1
        )

fig.write_image(f"recreated_boxplot_{{feature_id}}.svg")
# fig.show()
"""
    return script




if __name__ == "__main__":
    fbmn_task_id = "58e0e2959ec748049cb2c5f8bb8b87dc"
    cmmc_task_id = "21c17a8de65041369d607493140a367f"
    metadata_file = './data/metadata_quinn2020.tsv'

    enriched_df = fetch_enriched_results(cmmc_task_id)

    metadata_df = pd.read_csv(metadata_file, sep='\t')
    quant_df = fbmn_quant_download_wrapper(fbmn_task_id)
