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


@dataclass
class FilterResult:
    data: pd.DataFrame
    filtered: bool
    filters: str


def render_filter_options(merged_df, first_option, second_option, key: str) -> FilterResult:
    filter_on = st.checkbox("Use column and value filters", key=key)
    filter_str = ''

    if filter_on:
        filter_by = "source" if "input_source" in first_option else "origin"

        bool_matrix_df = prepare_dataframe(
            merged_df,
            by=filter_by
        )
        target_set = second_option
        matches_df = find_exact_matches(bool_matrix_df, target_set)
        merged_df = merged_df.iloc[matches_df.index]

        # Build filter strings - adjust based on what you want to capture
        filter_str = f"filtered by {filter_by} - values: {','.join(target_set)}"

    return FilterResult(
        data=merged_df,
        filtered=filter_on,
        filters=filter_str
    )


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


def create_pdf_download_button(
        data_df,
        feat_id_dict,
        plot_function,
        plot_params,
        button_key,
        button_label="Download all plots",
        button_icon=":material/download:",
        file_prefix="boxplots_all_features"
):
    """
    Creates a PDF download button that generates plots for all feature IDs.

    Parameters:
    - data_df: The dataframe containing the data
    - feat_id_dict: Dictionary of feature IDs and their info
    - plot_function: The plotting function to call (e.g., box_plot.plot_boxplots_by_group)
    - plot_params: Dictionary of parameters to pass to the plotting function
    - button_label: Label for the download button
    - file_prefix: Prefix for the downloaded PDF filename

    Returns:
    - None (displays the download button in Streamlit)
    """

    if st.button(button_label, type="secondary", icon=button_icon, key=button_key):

        if len(feat_id_dict) == 0:
            st.warning("No feature IDs available for PDF generation")
            return

        try:
            # Create a temporary file for the PDF
            with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tmp_file:
                pdf_path = tmp_file.name

            # Generate PDF with all plots
            with PdfPages(pdf_path) as pdf:
                success_count = 0

                for feat_id, feat_info in feat_id_dict.items():
                    try:
                        # Update plot parameters with current feature ID
                        current_params = plot_params.copy()
                        current_params['feature_id'] = int(feat_id)
                        current_params['df_quant_merged'] = data_df

                        # Generate plotly figure
                        fig, _ = plot_function(**current_params)

                        # Convert plotly figure to matplotlib and save to PDF
                        img_bytes = pio.to_image(fig, format="png", width=1200, height=800)

                        # Create matplotlib figure
                        plt.figure(figsize=(12, 8))
                        plt.imshow(plt.imread(BytesIO(img_bytes)))
                        plt.axis('off')
                        plt.title(f"Feature ID {feat_id}: {feat_info}",
                                  fontsize=14, pad=20)

                        # Save to PDF
                        pdf.savefig(bbox_inches='tight', dpi=300)
                        plt.close()
                        success_count += 1

                    except Exception as e:
                        st.warning(f"Could not generate plot for Feature ID {feat_id}: {str(e)}")
                        continue

            if success_count > 0:
                # Read the PDF file and provide download
                with open(pdf_path, "rb") as pdf_file:
                    pdf_bytes = pdf_file.read()

                # Create filename based on plot parameters
                filename_suffix = plot_params.get('filename_suffix', 'plots')
                filename = f"{file_prefix}_{filename_suffix}.pdf"

                st.download_button(
                    label="📁 Download PDF",
                    data=pdf_bytes,
                    file_name=filename,
                    mime="application/pdf"
                )

                st.success(f"PDF generated with {success_count} plots!")
            else:
                st.error("No plots could be generated. Please check your selections.")

            # Clean up temporary file
            os.unlink(pdf_path)

        except Exception as e:
            st.error(f"Error generating PDF: {str(e)}")


# Usage for the first section (main box plots with groups1 and groups2):
def add_pdf_download_boxplots(merged_data, feat_id_dict, groups2, groups1, column1, column2,
                              boxp_filter_string, color_mapping):
    """For the first section with group 1 and group 2 parameters"""

    plot_params = {
        'groups2': [groups2] if groups2 else None,
        'groups1': groups1,
        'column1': column1,
        'column2': column2,
        'informations': boxp_filter_string,
        'color_mapping': color_mapping,
    }

    create_pdf_download_button(
        data_df=merged_data,
        feat_id_dict=feat_id_dict,
        plot_function=box_plot.plot_boxplots_by_group,
        plot_params=plot_params,
        button_key="boxplot_download",
        file_prefix="boxplots_analysis"
    )


# Usage for the second section (overview with single group):
def add_pdf_download_overview(data_overview_df, feat_id_dict, group_by, column_select, filter_string, color_mapping):
    """For the second section with single group comparison"""

    plot_params = {
        'groups1': group_by,
        'column1': column_select,
        'informations': filter_string,
        'color_mapping': color_mapping,
    }

    create_pdf_download_button(
        data_df=data_overview_df,
        feat_id_dict=feat_id_dict,
        plot_function=box_plot.plot_boxplots_by_group,
        plot_params=plot_params,
        button_key="overview_donwload",
        file_prefix="boxplots_overview"
    )


def load_uploaded_file_df(uploaded_file):
    if uploaded_file.name.endswith(".csv"):
        loaded_file_df = pd.read_csv(uploaded_file)
    elif uploaded_file.name.endswith(".tsv"):
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
    # First, identify the data rows we need (those containing mzML/mzXML)
    data_mask = df_quant.columns.str.contains("mzML|mzXML", case=False, na=False)
    if not data_mask.any():
        # If no columns match, check the index after transpose
        t_df_quant = df_quant.set_index("row ID").T
        data_mask = t_df_quant.index.str.contains("mzML|mzXML", case=False, na=False)
        t_df_quant = t_df_quant[data_mask].reset_index()
        t_df_quant.columns.name = None
        t_df_quant = t_df_quant.rename(columns={"index": "filename"})
    else:
        # Standard transpose operation
        t_df_quant = df_quant.set_index("row ID").T
        data_mask = t_df_quant.index.str.contains("mzML|mzXML", case=False, na=False)
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
    filename_pattern = r"\.mzML|\.mzXML"
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
        st.markdown(f"- [Contribute depositing more information for {feature_id.split(':')[1]} on CMMC-kb]({url})")
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
        inchikey = smiles_to_inchikey(smiles)
        enrichment_date_ = st.session_state['enrichment_date']
        smiles_svg = smiles_to_svg(smiles, (500, 500))
        if smiles_svg:
            st.image(smiles_svg)
            st.markdown(f"""Data as of {enrichment_date_}""", unsafe_allow_html=True)
            st.markdown(f"[View latest data](https://cmmc-kb.gnps2.org/structurepage/?inchikey={inchikey})")
        else:
            st.info("No valid SMILES string available to generate structure image.")
            st.markdown(f"""Data as of {enrichment_date_}""", unsafe_allow_html=True)
        st.markdown("<br>".join(text_info), unsafe_allow_html=True)
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


def params_to_binary_string(params):
    """Convert parameters to a compressed base64 string"""
    # Convert to JSON string
    json_str = json.dumps(params, separators=(',', ':'))  # Compact JSON

    # Compress the JSON string
    buffer = BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode='wb') as f:
        f.write(json_str.encode('utf-8'))

    # Encode as base64
    compressed_data = buffer.getvalue()
    b64_string = base64.b64encode(compressed_data).decode('ascii')

    return b64_string


def binary_string_to_params(b64_string):
    """Convert base64 string back to parameters"""
    try:
        # Decode base64
        compressed_data = base64.b64decode(b64_string.encode('ascii'))

        # Decompress
        buffer = BytesIO(compressed_data)
        with gzip.GzipFile(fileobj=buffer, mode='rb') as f:
            json_str = f.read().decode('utf-8')

        # Parse JSON
        params = json.loads(json_str)
        return params, None

    except Exception as e:
        return None, str(e)


if __name__ == "__main__":
    fbmn_task_id = "58e0e2959ec748049cb2c5f8bb8b87dc"
    cmmc_task_id = "21c17a8de65041369d607493140a367f"
    metadata_file = './data/metadata_quinn2020.tsv'

    enriched_df = fetch_enriched_results(cmmc_task_id)

    metadata_df = pd.read_csv(metadata_file, sep='\t')
    quant_df = fbmn_quant_download_wrapper(fbmn_task_id)
