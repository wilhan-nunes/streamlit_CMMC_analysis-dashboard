import requests
import pandas as pd
from typing import Literal


def fetch_enriched_results(task_id: str) -> pd.DataFrame:
    """
    Download enriched results as a tsv file from GNPS2 Enrichment workflow.

    :param task_id: GNPS2 Enrichment workflow task ID
    :return: pd.DataFrame containing enriched results
    """
    url = f"https://gnps2.org/resultfile?task={task_id}&file=nf_output/cmmc_results/cmmc_enriched_results.tsv"
    df = pd.read_csv(url, sep='\t', low_memory=False)
    return df


def fetch_phylogeny_results(task_id: str) -> pd.DataFrame:
    """
    Fetch phylogeny results from GNPS2.

    :param task_id: GNPS2 Enrichment workflow task ID
    :return: pd.DataFrame containing phylogeny results
    """
    url = f"https://gnps2.org/resultfile?task={task_id}&file=nf_output/cmmc_results/cmmc_taxonomy.tsv"
    df = pd.read_csv(url, sep='\t', low_memory=False)
    return df


def fetch_cmmc_graphml(task_id: str, graphml_path='data/network.graphml'):
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
            with open(filepath, 'wb') as f:
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


if __name__ == '__main__':
    # AGP example
    cmmc_task_id = '1715c16a223e47c98d0a70a26ef6f8ef'
    fbmn_task_id = '0a5bc6c69d7c4827824c6329804f2c12'
    # http://localhost:8501/?cmmc_task_id=1715c16a223e47c98d0a70a26ef6f8ef&fbmn_task_id=0a5bc6c69d7c4827824c6329804f2c12

    # # 3D Mice
    # fbmn_task_id = '58e0e2959ec748049cb2c5f8bb8b87dc'
    # cmmc_task_id = '21c17a8de65041369d607493140a367f'
    # metadata_file = 'data/metadata_quinn2020.txt'

    enriched_result = fetch_enriched_results(cmmc_task_id)
    print(f"Enriched Results for Task ID {cmmc_task_id}:\n", enriched_result.head())

    phylogeny_result = fetch_phylogeny_results(cmmc_task_id)
    print(f"Phylogeny Results for Task ID {cmmc_task_id}:\n", phylogeny_result.head())

    graphml_file_name = fetch_cmmc_graphml(cmmc_task_id, graphml_path=f'data/{cmmc_task_id}_network.graphml')
    print(f"GraphML file saved at: {graphml_file_name}")

    fbmn_quant_table = fetch_file(fbmn_task_id, 'mouse_fbmn_quant_table.csv', type='quant_table')
    print(f"FBMN Quantification Table saved at: {fbmn_quant_table}")
    quant_table = pd.read_csv(fbmn_quant_table)
    # quant_table = pd.read_csv('./data/American_Gut_Project_-_all_lib-0a5bc6c69d7c4827824c6329804f2c12-featuretable_reformated.csv')

    # metadata AGP
    metadata = fetch_file(fbmn_task_id, 'fbmn_metadata.tsv', type='annotation_table')
    metadata_df = pd.read_csv(metadata, sep='\t')

    # #metadata 3D Mouse
    # metadata_file = 'data/metadata_quinn2020.txt'
    # metadata_df = pd.read_csv(metadata_file, sep='\t')

