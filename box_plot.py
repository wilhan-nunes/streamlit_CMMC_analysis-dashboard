import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from utils import fetch_enriched_results


def prepare_lcms_data(
        df_quant: pd.DataFrame, df_metadata: pd.DataFrame, cmmc_results: pd.DataFrame
):
    df_quant = df_quant.copy()
    microbial_scans = cmmc_results["query_scan"].tolist()

    df_quant = df_quant[df_quant["row ID"].isin(microbial_scans)]

    # Remove "Peak area" from column1 names
    df_quant.columns = df_quant.columns.str.replace(" Peak area", "")
    # Transpose dataframe
    t_df_quant = df_quant.transpose()
    # Set column1 names from row.ID
    t_df_quant.columns = df_quant["row ID"]
    # Remove first 3 rows
    t_df_quant = t_df_quant[t_df_quant.index.str.contains("mzML|mzXML", case=False, na=False)]
    # Reset index to create filename column1
    t_df_quant = t_df_quant.reset_index()
    t_df_quant = t_df_quant.rename(columns={"index": "filename"})
    # Melt the dataframe to create long format
    df_quant_long = pd.melt(
        t_df_quant, id_vars=["filename"], var_name="featureID", value_name="Abundance"
    )
    # Merge dataframes
    df_quant_merged = pd.merge(df_metadata, df_quant_long, on="filename")
    # merge with cmmc results
    cmmc_results = cmmc_results.rename(columns={"query_scan": "featureID"})
    df_quant_merged = pd.merge(
        df_quant_merged, cmmc_results, on="featureID", how="left"
    )
    return df_quant_merged


def plot_boxplots_by_group(
        df_quant_merged: pd.DataFrame,
        groups1: list = None,
        groups2: list = None,
        feature_id: int = None,
        column1: str = None,
        column2: str = None,
):
    """
    Plots boxplots of 'Abundance' for selected groups in "column" using Plotly Express.
    Args:
        df_quant_merged (pd.DataFrame): The merged dataframe.
        groups1 (list): List of group names to filter and plot.
        feature_id (int): The feature ID to filter by.
        column1 (str): The column name to group by.
    Returns:
        plotly.graph_objects.Figure: The figure object for further use (e.g., in Streamlit).
    """

    # Optionally prefilter by groups1 and column1 if provided
    filtered_df = df_quant_merged.copy()
    if groups1 and column1:
        filtered_df = filtered_df[filtered_df[column1].isin(groups1)]
    if groups2 and column2:
        filtered_df = filtered_df[filtered_df[column2].isin(groups2)]
    if feature_id:
        filtered_df = filtered_df[filtered_df["featureID"] == feature_id]

    # Check if filtered data is empty
    if filtered_df.empty:
        print(f"Warning: No data found for feature_id {feature_id} and groups1 {groups1}")
        # Create empty plot
        fig = go.Figure()
        fig.add_annotation(
            text="No data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, xanchor='center', yanchor='middle',
            showarrow=False, font=dict(size=16)
        )
        fig.update_layout(
            title=f"No data for Feature ID: {feature_id}",
            xaxis_title=column1,
            yaxis_title="Abundance",
            width=800, height=500
        )
        return fig

    # Get the title
    title_data = df_quant_merged[df_quant_merged["featureID"] == feature_id]["input_name"]
    if not title_data.empty:
        title = title_data.iloc[0]
    else:
        title = f"Feature ID: {feature_id}"

    try:
        # Create boxplot with Plotly Express
        fig = px.box(
            filtered_df,
            x=column1,
            y="Abundance",
            category_orders={column1: groups1},  # This maintains the order
            title=str(title).capitalize(),
            width=800,
            height=500,
            points="all",
        )

        # Remove fliers (outliers)
        fig.update_traces(boxpoints='all', marker=dict(opacity=0.6))

        # Update layout for better appearance
        fig.update_layout(
            xaxis_title=column1,
            yaxis_title="Abundance",
            showlegend=False
        )

    except Exception as e:
        print(f"Error creating boxplot: {e}")
        # Create error plot
        fig = go.Figure()
        fig.add_annotation(
            text=f"Error: {str(e)}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, xanchor='center', yanchor='middle',
            showarrow=False, font=dict(size=16, color="red")
        )
        fig.update_layout(
            title=f"Error plotting Feature ID: {feature_id}",
            xaxis_title=column1,
            yaxis_title="Abundance",
            width=800, height=500
        )

    return fig

if __name__ == "__main__":
    # Example usage
    from utils import fetch_file

    fbmn_task_id = "58e0e2959ec748049cb2c5f8bb8b87dc"
    cmmc_task_id = "21c17a8de65041369d607493140a367f"

    cmmc_result = fetch_enriched_results(cmmc_task_id)
    quant_file = fetch_file(fbmn_task_id, "mouse_metadata.csv", "quant_table")

    df_quant = pd.read_csv(quant_file)
    metadata = pd.read_csv("data/metadata_quinn2020.tsv", sep="\t")

    merged_df = prepare_lcms_data(df_quant, metadata, cmmc_result)

    # This will be also the order of the boxplots
    keywords = ["Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Stool"]
    id_num = 29384
    #
    boxplot_fig = plot_boxplots_by_group(merged_df, keywords, id_num)
    boxplot_fig.savefig(f"boxplot_by_body_part_{id_num}.png")
