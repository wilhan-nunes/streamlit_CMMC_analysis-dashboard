import pandas as pd
from upsetplot import UpSet, from_indicators
import matplotlib.pyplot as plt
from typing import Literal


def generate_upset_plot(enrich_df: pd.DataFrame, by: Literal["source", "origin"]) -> plt.Figure:
    """
    Generate an UpSet plot from the input source or origin data in the enrichment DataFrame.
    :param enrich_df:
    :param by: "origin" or "source" to determine which column to use for the UpSet plot.
    :return: figure object containing the UpSet plot.
    """
    # Clean and standardize the 'input_source' column
    enrich_df['input_source_clean'] = (
        enrich_df['input_source']
        .fillna('')
        .str.replace(r'\s+and\s+', ';', regex=True)  # Normalize 'and' to ';'
        .str.split(';')
    )
    # Extract all unique source categories
    all_sources = sorted({source.strip() for sources in enrich_df['input_source_clean'] for source in sources if source})
    # Create boolean indicator columns for each source
    df_indicators = pd.DataFrame()
    for source in all_sources:
        df_indicators[source] = enrich_df['input_source_clean'].apply(lambda x: source in x)
    # Omit rows where all indicators are False
    df_indicators = df_indicators[df_indicators.any(axis=1)]
    # Generate data for the UpSet plot
    upset_data = from_indicators(df_indicators.columns.tolist(), df_indicators)
    # Plot the UpSet diagram on a specific figure and axes
    fig, ax = plt.subplots(figsize=(10, 6))
    UpSet(upset_data, subset_size='count').plot(fig)
    fig.suptitle("UpSet Plot of input_source Categories", y=1.05)
    plt.tight_layout()
    return fig
