import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from utils import *


def plot_boxplots_by_group(
        df_quant_merged: pd.DataFrame,
        groups1: list = None,
        groups2: list = None,
        feature_id: int = None,
        column1: str = None,
        column2: str = None,
        informations: str = None,
        color_mapping: dict = None):
    """
    Plots boxplots of 'Peak Area' for selected groups in "column" using Plotly Express.
    Args:
        df_quant_merged (pd.DataFrame): The merged dataframe.
        groups1 (list): List of group names to filter and plot.
        feature_id (int): The feature ID to filter by.
        column1 (str): The column name to group by.
        color_mapping (dict): Optional mapping for colors.
    Returns:
        plotly.graph_objects.Figure: The figure object for further use (e.g., in Streamlit).
        pd.DataFrame: filtered df used to generate the boxplot
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
            yaxis_title="Peak Area",
            width=800, height=500
        )
        return fig, filtered_df

    # Get the title
    title_data = df_quant_merged[df_quant_merged["featureID"] == feature_id]["input_name"]
    if not title_data.empty:
        title = str(title_data.iloc[0])
    else:
        title = f"Feature ID: {feature_id}"
    if informations:
        title = title + " | " + informations

    try:
        # Create boxplot with Plotly Express
        fig = px.box(
            filtered_df,
            x=column1,
            y="Peak Area",
            category_orders={column1: groups1},  # This maintains the order
            color=column1,
            color_discrete_sequence=px.colors.qualitative.Set1,
            color_discrete_map=None if not color_mapping else color_mapping,
            width=800,
            height=500,
            points="all",
            hover_name='filename'
        )

        # Remove fliers (outliers)
        fig.update_traces(boxpoints='all',
                          marker=dict(opacity=0.6),
                          pointpos=0,
                          jitter=0.3
                          )

        # Update layout for better appearance
        fig.update_layout(
            title=str(title).capitalize(),
            xaxis_title=column1,
            yaxis_title="Peak Area",
            showlegend=False,
            yaxis=dict(exponentformat="power", showexponent="all"),
        )

        group_counts = filtered_df.groupby(column1)['Peak Area'].count()
        tickvals = [g for g in groups1 if g in group_counts.index]
        ticktext = [f"{g}<br>(N={group_counts[g]})" for g in tickvals]

        fig.update_xaxes(tickvals=tickvals, ticktext=ticktext)

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
            yaxis_title="Peak Area",
            width=800, height=500
        )

    return fig, filtered_df


def insert_plot_download_buttons(plot_data_df, identifier, svg_data, key_prefix):
    col1, col2, _, _ = st.columns(4)
    col1.download_button(
        label=":material/download: Download Plot as SVG",
        data=svg_data,
        file_name=f"network_{identifier}.svg",
        mime="image/svg+xml",  # Set the MIME type to SVG
        key=f'{key_prefix}_svg_download'
    )
    col2.download_button(
        label=":material/download: Download Plot Data",
        data=plot_data_df.to_csv(index=False, sep='\t'),
        file_name=f"box_plot_data_{identifier}.tsv",
        mime="text/tsv",
        key=f'{key_prefix}_tsv_download'
    )


if __name__ == "__main__":
    from utils import fetch_enriched_results, prepare_lcms_data

    # fbmn_task_id = "58e0e2959ec748049cb2c5f8bb8b87dc"
    # cmmc_task_id = "21c17a8de65041369d607493140a367f"

    fbmn_task_id = "0f32c284e015471c96e39bbbefe8f2d2"
    cmmc_task_id = "754604934a5a4730a255abb0c96c7111"

    cmmc_result = fetch_enriched_results(cmmc_task_id)
    df_quant = fbmn_quant_download_wrapper(fbmn_task_id)

    # metadata = pd.read_csv("data/metadata_quinn2020.tsv", sep="\t")
    metadata = pd.read_csv("data/merged_metadata.tsv", sep="\t")

    merged_df = prepare_lcms_data(df_quant, metadata, cmmc_result, include_all_scans=False)

    keywords = ["Stomach", "Duodenum", "Jejunum", "Ileum", "Cecum", "Colon", "Stool"]
    id_num = 29384

    boxplot_fig, box_plot_data = plot_boxplots_by_group(
        merged_df,
        groups1=keywords,
        feature_id=id_num,
        column1="ATTRIBUTE_UBERONBodyPartName"  # Replace with the actual column name in your metadata
    )
    boxplot_fig.write_image(f"boxplot_by_body_part_{id_num}.png")
