import io
import json
import os
import tempfile
import urllib.parse
from dataclasses import dataclass
from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd
import plotly.io as pio
import requests
import streamlit as st
from gnpsdata import taskinfo, taskresult, workflow_fbmn
from matplotlib.backends.backend_pdf import PdfPages
from rdkit import Chem
from rdkit.Chem import Draw

import box_plot


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
def fetch_cmmc_graphml(task_id: str, graphml_path="data/network.graphml"):
    """
    Fetch CMMC graphml results from GNPS2.

    :param task_id: GNPS2 Enrichment workflow task ID
    :param graphml_path: Path to save the graphml file
    :return: path to CMMC graphml results
    """
    url = taskresult.determine_gnps2_resultfile_url(task_id,"nf_output/gnps_network/network.graphml")
    with requests.get(url) as response:
        if response.status_code == 200:
            filepath = graphml_path
            with open(filepath, "wb") as f:
                f.write(response.content)
            return filepath
        else:
            raise Exception(f"Failed to fetch graphml file: {response.status_code}")


# this is used for the FBMN files
def fetch_file(
        task_id: str, type: Literal["quant_table", "annotation_table"]
) -> pd.DataFrame:
    """
    Fetches a file from a given task ID and loads it into a pandas DataFrame.

    :param task_id: The task ID to construct the file URL.
    :param file_name: The name of the file to fetch. Must be one of the predefined options.
    :param type: The type of file to fetch. Must be one of "quant_table" or "library_search_table".
    :returns: The path to the downloaded file.
    """

    df = None
    if type == "annotation_table":
        df = workflow_fbmn.get_metadata_dataframe(task_id, gnps2=True)
    elif type == "quant_table":
        df = workflow_fbmn.get_quantification_dataframe(task_id, gnps2=True)

    return df


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
                        plt.imshow(plt.imread(io.BytesIO(img_bytes)))
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
            if validation_str not in workflow_name:
                st.error(f"The provided task ID is for {workflow_name}", icon=":material/error:")
        except Exception:
            st.error(f"Failed to fetch task information. Is this a valid task id?", icon=":material/dangerous:")
    elif task_id:
        st.warning("Task ID must be exactly 32 characters long.")


if __name__ == "__main__":
    fbmn_task_id = "58e0e2959ec748049cb2c5f8bb8b87dc"
    cmmc_task_id = "21c17a8de65041369d607493140a367f"
    metadata_file = './data/metadata_quinn2020.tsv'

    enriched_df = fetch_enriched_results(cmmc_task_id)

    graphml_file_name = fetch_cmmc_graphml(
        cmmc_task_id, graphml_path=f"data/{cmmc_task_id}_network.graphml"
    )

    metadata_df = pd.read_csv(metadata_file, sep='\t')
    quant_df = fetch_file(fbmn_task_id, type="quant_table")
