import warnings
from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet, from_indicators

# Suppress all warnings from upsetplot
warnings.filterwarnings("ignore", module="upsetplot")


def generate_upset_plot(
    enrich_df: pd.DataFrame, by: Literal["source", "origin"]
) -> plt.Figure:
    """
    Generate an UpSet plot from the input source or origin data in the enrichment DataFrame.

    :param enrich_df: DataFrame containing enrichment data with either 'input_source' or 'input_molecule_origin'
    :param by: "origin" or "source" to determine which column1 to use for the UpSet plot.
    :return: figure object containing the UpSet plot.
    """
    # Select the column1 to process
    print(f"Generating UpSet plot by '{by}'...")
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

    # Extract all unique categories
    all_categories = sorted(
        {item.strip() for items in enrich_df["input_clean"] for item in items if item}
    )

    # Create boolean indicator columns for each category
    df_indicators = pd.DataFrame()
    for category in all_categories:
        df_indicators[category] = enrich_df["input_clean"].apply(
            lambda x: category in x
        )

    # Omit rows where all indicators are False
    df_indicators = df_indicators[df_indicators.any(axis=1)]

    # Generate data for the UpSet plot
    upset_data = from_indicators(df_indicators.columns.tolist(), df_indicators)

    # Plot the UpSet diagram
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_axis_off()
    UpSet(
        upset_data, subset_size="count", sort_by="cardinality", show_counts=True
    ).plot(fig)
    fig.suptitle(f"UpSet Plot for input_{by} Categories", y=1.05)
    for ax_ in fig.axes:
        ax_.grid(axis="x", visible=False)

    # Convert the figure to SVG and return as a string
    import io

    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    svg = buf.getvalue()
    buf.close()
    plt.close(fig)
    return svg
