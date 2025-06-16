import io
import os
import tempfile
from dataclasses import dataclass
from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd
import plotly.io as pio
import requests
import streamlit as st
from matplotlib.backends.backend_pdf import PdfPages

import box_plot


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
        filter_str = f"filtered by {filter_by} - values: {",".join(target_set)}"

    return FilterResult(
        data=merged_df,
        filtered=filter_on,
        filters=filter_str
    )


def fetch_enriched_results(task_id: str) -> pd.DataFrame:
    """
    Download enriched results as a tsv file from GNPS2 Enrichment workflow.

    :param task_id: GNPS2 Enrichment workflow task ID
    :return: pd.DataFrame containing enriched results
    """
    url = f"https://gnps2.org/resultfile?task={task_id}&file=nf_output/cmmc_results/cmmc_enriched_results.tsv"
    df = pd.read_csv(url, sep="\t", low_memory=False)
    return df


def fetch_phylogeny_results(task_id: str) -> pd.DataFrame:
    """
    Fetch phylogeny results from GNPS2.

    :param task_id: GNPS2 Enrichment workflow task ID
    :return: pd.DataFrame containing phylogeny results
    """
    url = f"https://gnps2.org/resultfile?task={task_id}&file=nf_output/cmmc_results/cmmc_taxonomy.tsv"
    df = pd.read_csv(url, sep="\t", low_memory=False)
    return df


def fetch_cmmc_graphml(task_id: str, graphml_path="data/network.graphml"):
    """
    Fetch CMMC graphml results from GNPS2.

    :param task_id: GNPS2 Enrichment workflow task ID
    :param graphml_path: Path to save the graphml file
    :return: path to CMMC graphml results
    """
    url = f"https://gnps2.org/resultfile?task={task_id}&file=nf_output/gnps_network/network.graphml"
    with requests.get(url) as response:
        if response.status_code == 200:
            filepath = graphml_path
            with open(filepath, "wb") as f:
                f.write(response.content)
            return filepath
        else:
            raise Exception(f"Failed to fetch graphml file: {response.status_code}")


def fetch_file(
        task_id: str, file_name: str, type: Literal["quant_table", "annotation_table"]
) -> str:
    """
    Fetches a file from a given task ID and loads it into a pandas DataFrame.

    :param task_id: The task ID to construct the file URL.
    :param file_name: The name of the file to fetch. Must be one of the predefined options.
    :param type: The type of file to fetch. Must be one of "quant_table" or "library_search_table".
    :returns: The path to the downloaded file.
    """
    if type == "annotation_table":
        input_url = f"https://gnps2.org/resultfile?task={task_id}&file=nf_output/library/merged_results_with_gnps.tsv"
    elif type == "quant_table":
        input_url = f"https://gnps2.org/result?task={task_id}&viewname=quantificationdownload&resultdisplay_type=task"
    response = requests.get(input_url)
    response.raise_for_status()  # Raise an error for failed requests
    output_file_path = f"data/{file_name}"

    with open(output_file_path, "w") as f:
        f.write(response.text)

    return output_file_path


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
        button_label="📥 Download All Plots as PDF",
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

    if st.button(button_label, type="secondary", key=button_key):

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
                        fig = plot_function(**current_params)

                        # Convert plotly figure to matplotlib and save to PDF
                        img_bytes = pio.to_image(fig, format="png", width=1200, height=800)

                        # Create matplotlib figure
                        plt.figure(figsize=(12, 8))
                        plt.imshow(plt.imread(io.BytesIO(img_bytes)))
                        plt.axis('off')
                        plt.title(f"Feature ID {feat_id}: {feat_info.get('input_name', 'Unknown')}",
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
                              boxp_filter_string):
    """For the first section with group 1 and group 2 parameters"""

    plot_params = {
        'groups2': [groups2] if groups2 else None,
        'groups1': groups1,
        'column1': column1,
        'column2': column2,
        'informations': boxp_filter_string,
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
def add_pdf_download_overview(data_overview_df, feat_id_dict, group_by, column_select, filter_string):
    """For the second section with single group comparison"""

    plot_params = {
        'groups1': group_by,
        'column1': column_select,
        'informations': filter_string,
    }

    create_pdf_download_button(
        data_df=data_overview_df,
        feat_id_dict=feat_id_dict,
        plot_function=box_plot.plot_boxplots_by_group,
        plot_params=plot_params,
        button_key="overview_donwload",
        file_prefix="boxplots_overview"
    )


if __name__ == "__main__":
    # AGP example
    # cmmc_task_id = '1715c16a223e47c98d0a70a26ef6f8ef'
    # fbmn_task_id = '0a5bc6c69d7c4827824c6329804f2c12'
    # http://localhost:8501/?cmmc_task_id=1715c16a223e47c98d0a70a26ef6f8ef&fbmn_task_id=0a5bc6c69d7c4827824c6329804f2c12

    # # 3D Mice
    fbmn_task_id = "58e0e2959ec748049cb2c5f8bb8b87dc"
    cmmc_task_id = "21c17a8de65041369d607493140a367f"
    metadata_file = 'data/metadata_quinn2020.txt'

    enriched_df = fetch_enriched_results(cmmc_task_id)

    # phylogeny_result = fetch_phylogeny_results(cmmc_task_id)

    graphml_file_name = fetch_cmmc_graphml(
        cmmc_task_id, graphml_path=f"data/{cmmc_task_id}_network.graphml"
    )

    # fbmn_quant_table = fetch_file(
    #     fbmn_task_id, "mouse_fbmn_quant_table.csv", type="quant_table"
    # )
    # print(f"FBMN Quantification Table saved at: {fbmn_quant_table}")
    # quant_table = pd.read_csv(fbmn_quant_table)

    # #metadata 3D Mouse
    metadata_file = 'data/metadata_quinn2020.txt'
    metadata_df = pd.read_csv(metadata_file, sep='\t')
