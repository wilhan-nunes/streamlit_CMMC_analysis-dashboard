import networkx as nx
import plotly.graph_objects as go


def plot_cluster_by_node(
        G,
        node_id,
        annotate_nodes,
        nodes_info="id",
        width=1000,
        height=700,
        layout="kamada",
        node_colors_dict=None,
        show_delta_annotation=False,
):
    """
    Find the component ID for a given node and plot that entire cluster.

    Parameters:
    G: NetworkX graph
    node_id: The node ID to look up
    annotate_nodes: list with all nodes to color in the network
    width, height: plot dimensions

    Returns:
    fig: Plotly figure object
    """

    # Check if node exists
    if node_id not in G.nodes():
        return None

    if node_colors_dict is None:
        node_colors_dict = {
            "queried_node": "rgba(210,55,44,1)",
            "cmmc_match": "rgba(44, 146, 31, 1)",
            "fbmn_match": "rgba(199, 133, 7, 1)",
            "unannotated": "rgba(75,125,180,1)",
        }

    # Get the component ID for this node
    component_id = G.nodes[node_id].get("component")

    if component_id is None:
        return None

    # Find all nodes in this component
    cluster_nodes = [
        node for node in G.nodes() if G.nodes[node].get("component") == component_id
    ]

    cluster_mz_values = [G.nodes[node].get("mz", None) for node in cluster_nodes]

    # Find all edges in this component
    cluster_edges = [
        edge
        for edge in G.edges()
        if edge[0] in cluster_nodes and edge[1] in cluster_nodes
    ]

    # Create subgraph for this component
    subgraph = G.subgraph(cluster_nodes)

    # Choose layout
    if layout == "spring":
        pos = nx.spring_layout(subgraph, k=2, iterations=100, seed=42)
    else:
        pos = nx.kamada_kawai_layout(subgraph)

    # Prepare node coordinates
    x_nodes = [pos[node][0] for node in cluster_nodes]
    y_nodes = [pos[node][1] for node in cluster_nodes]

    # Prepare edge coordinates for lines (visual edges)
    x_edges_line, y_edges_line = [], []
    for edge in subgraph.edges():
        x_edges_line.extend([pos[edge[0]][0], pos[edge[1]][0], None])
        y_edges_line.extend([pos[edge[0]][1], pos[edge[1]][1], None])

    # Prepare invisible points along edges for hover
    x_edges_hover, y_edges_hover = [], []
    edge_hover_text = []
    annotation_edge_delta = []
    annotation_edge_delta_x = []
    annotation_edge_delta_y = []

    for edge in cluster_edges:
        edge_attr = G.edges[edge]

        # Create edge hover text
        edge_text = f"<b>Edge: {edge[0]} → {edge[1]}</b><br>"
        edge_text += f"Δm/z: {edge_attr.get('deltamz', 'N/A')}<br>"

        # Add any other edge attributes you want to display
        if "score" in edge_attr:
            edge_text += f"Cosine Score: {edge_attr['score']:.3f}<br>"
        if "mass_diff" in edge_attr:
            edge_text += f"Mass Diff: {edge_attr['mass_diff']:.4f}<br>"

        # Get positions of the two nodes
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]

        # Create multiple points along the edge (10 points for good coverage)
        num_points = 10
        for i in range(num_points):
            t = i / (num_points - 1)  # Parameter from 0 to 1
            x_point = x0 + t * (x1 - x0)
            y_point = y0 + t * (y1 - y0)

            x_edges_hover.append(x_point)
            y_edges_hover.append(y_point)
            edge_hover_text.append(edge_text)

            if i == 5:
                # Add a detailed annotation for the middle point of the edge
                annotation_edge_delta.append(
                    f"Edge: {edge[0]} → {edge[1]}<br>"
                    f"Δm/z: {edge_attr.get('deltamz', 'N/A')}<br>"
                )
                annotation_edge_delta_x.append(x_point)
                annotation_edge_delta_y.append(y_point)
    # Get node data for hover text
    hover_text = []
    node_colors = []

    for node in cluster_nodes:
        attr = G.nodes[node]

        # Create hover text
        text = f"<b>Node {node}</b><br>"
        text += f"m/z: {attr.get('mz', 'N/A')}<br>"
        text += f"RT: {attr.get('rt', 'N/A')} min<br>"

        if "library_compound_name" in attr and attr["library_compound_name"]:
            text += f"<b>Compound:</b><br>{attr['library_compound_name']}<br>"
        else:
            text += f"<i>Unidentified</i><br>"

        hover_text.append(text)

        # Color the queried node differently
        if node == node_id:
            node_colors.append(
                node_colors_dict.get("queried_node", "rgba(210,55,44,1)")
            )  # Highlight the queried node
        elif node in annotate_nodes:
            node_colors.append(
                node_colors_dict.get("cmmc_match", "rgba(44, 146, 31, 1)")
            )
        elif "library_compound_name" in attr:
            node_colors.append(
                node_colors_dict.get("fbmn_match", "rgba(199, 133, 7, 1)")
            )
        else:
            node_colors.append(
                node_colors_dict.get("unannotated", "rgba(75,125,180,1)")
            )

    # Create line trace for visual edges (no hover)
    line_trace = go.Scatter(
        x=x_edges_line,
        y=y_edges_line,
        line=dict(width=2, color="gray"),
        hoverinfo="skip",  # Skip hover for the line trace
        mode="lines",
        showlegend=False,
        name="Edge Lines",
    )

    # Create invisible hover points along edges
    edge_hover_trace = go.Scatter(
        x=x_edges_hover,
        y=y_edges_hover,
        mode="markers",
        marker=dict(
            size=8,  # Size of invisible markers
            color="rgba(0,0,0,0)",  # Completely transparent
            line=dict(width=0),
        ),
        hoverinfo="text",
        hovertext=edge_hover_text,
        showlegend=False,
        name="Edge Hover Points",
    )

    # Create node trace

    nodes_text = (
        [str(node) for node in cluster_nodes]
        if nodes_info == "id"
        else [str(node_mz) for node_mz in cluster_mz_values]
    )
    node_trace = go.Scatter(
        x=x_nodes,
        y=y_nodes,
        mode="markers+text",
        hoverinfo="text",
        hovertext=hover_text,
        marker=dict(size=30, color=node_colors, line=dict(width=2, color="white")),
        text=nodes_text,
        textposition="middle center",
        textfont=dict(size=18, color="black", shadow="0px -0px 2px white"),
        showlegend=False,
    )

    # Create delta annotation trace if requested
    delta_annotation_trace = None
    if show_delta_annotation:
        delta_annotation_trace = go.Scatter(
            x=annotation_edge_delta_x,
            y=annotation_edge_delta_y,
            mode="markers+text",
            marker=dict(size=0.1, color="rgba(0,0,0,0)"),  # Invisible markers
            hoverinfo="text",
            hovertext=annotation_edge_delta,
            text=annotation_edge_delta,
            textposition="top center",
            textfont=dict(size=14, color="black", shadow="0px -0px 2px white"),
            showlegend=False,
        )

    # Create figure with all traces
    all_traces = [line_trace, edge_hover_trace, node_trace]

    if delta_annotation_trace is not None:
        all_traces.append(delta_annotation_trace)

    fig = go.Figure(data=all_traces)

    fig.update_layout(
        # font_shadow='auto',
        margin=dict(l=20, r=20, t=20, b=80),
        showlegend=False,
        hovermode="closest",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor="white",
        width=width,
        height=height,
    )

    # Add cluster info
    identified = sum(
        1
        for node in cluster_nodes
        if "library_compound_name" in G.nodes[node]
        and G.nodes[node]["library_compound_name"]
    )

    info_text = f"Network {component_id}<br>"
    info_text += f"Nodes: {len(cluster_nodes)}<br>"
    info_text += f"Edges: {subgraph.number_of_edges()}<br>"
    info_text += f"Identified: {identified}<br>"
    info_text += f"<b>Queried node: {node_id}</b>"

    # # Add node color legend as annotation at the bottom
    legend_items = [
        (
            "Queried node",
            node_colors_dict.get("queried_node", "rgba(210,55,44,1)"),
        ),  # red
        (
            "CMMC Match",
            node_colors_dict.get("cmmc_match", "rgba(44, 146, 31, 1)"),
        ),  # green
        (
            "FBMN Match",
            node_colors_dict.get("fbmn_match", "rgba(199, 133, 7, 1)"),
        ),  # yellow
        (
            "Unidentified node",
            node_colors_dict.get("unannotated", "rgba(75,125,180,1)"),
        ),  # blue
    ]
    legend_html = "    ".join(
        f"<span style='color:{color}; font-weight:bold; font-size:24px;'>●</span> {label}"
        for label, color in legend_items
    )

    fig.add_annotation(
        x=0.5,
        y=-0.08,
        xref="paper",
        yref="paper",
        text=legend_html,
        showarrow=False,
        align="center",
        bgcolor="rgba(255,255,255,0.95)",
        borderwidth=1,
        font=dict(size=18),
        xanchor="center",
        yanchor="bottom",
    )

    return fig, info_text


def get_cluster_id(G, node_id):
    """
    Simple function to just get the component ID for a node.

    Parameters:
    G: NetworkX graph
    node_id: The node ID to look up (same as Feature ID)

    Returns:
    component_id: The component ID, or None if not found
    """
    if node_id not in G.nodes():
        return None

    component_id = G.nodes[node_id].get("component")

    return component_id


def quick_cluster_plot(G, node_id):
    """
    One-liner function: provide node ID, get the plot.
    """
    fig = plot_cluster_by_node(G, node_id)
    if fig:
        fig.show()
    return fig
