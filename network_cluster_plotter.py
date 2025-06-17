import networkx as nx
import plotly.graph_objects as go

def plot_cluster_by_node(G, node_id, width=1000, height=700, layout="kamada"):
    """
    Find the component ID for a given node and plot that entire cluster.
    
    Parameters:
    G: NetworkX graph
    node_id: The node ID to look up
    width, height: plot dimensions
    
    Returns:
    fig: Plotly figure object
    """
    
    # Check if node exists
    if node_id not in G.nodes():
        # print(f"❌ Node {node_id} not found in the network!")
        return None
    
    # Get the component ID for this node
    component_id = G.nodes[node_id].get('component')
    
    if component_id is None:
        # print(f"❌ Node {node_id} has no component information!")
        return None
    
    # print(f"🔍 Node {node_id} belongs to component {component_id}")
    
    # Find all nodes in this component
    cluster_nodes = [node for node in G.nodes() 
                    if G.nodes[node].get('component') == component_id]

    # Find all edges in this component
    cluster_edges = [edge for edge in G.edges()
                     if edge[0] in cluster_nodes and edge[1] in cluster_nodes]

    # print(f"📊 Component {component_id} contains {len(cluster_nodes)} nodes")
    
    # Create subgraph for this component
    subgraph = G.subgraph(cluster_nodes)
    
    # Choose layout
    if layout == 'spring':
        pos = nx.spring_layout(subgraph, k=2, iterations=100, seed=42)
    else:
        pos = nx.kamada_kawai_layout(subgraph)
    
    # Prepare node coordinates
    x_nodes = [pos[node][0] for node in cluster_nodes]
    y_nodes = [pos[node][1] for node in cluster_nodes]
    
    # Prepare edge coordinates
    x_edges, y_edges = [], []
    for edge in subgraph.edges():
        x_edges.extend([pos[edge[0]][0], pos[edge[1]][0], None])
        y_edges.extend([pos[edge[0]][1], pos[edge[1]][1], None])
    
    # Get node data for hover text
    hover_text = []
    node_colors = []
    
    for node in cluster_nodes:
        attr = G.nodes[node]
        
        # Create hover text
        text = f"<b>Node {node}</b><br>"
        text += f"m/z: {attr.get('mz', 'N/A')}<br>"
        text += f"RT: {attr.get('rt', 'N/A')} min<br>"
        
        if 'library_compound_name' in attr and attr['library_compound_name']:
            text += f"<b>Compound:</b><br>{attr['library_compound_name']}<br>"
        else:
            text += f"<i>Unidentified</i><br>"
            
        hover_text.append(text)
        
        # Color the queried node differently
        if node == node_id:
            node_colors.append('rgba(210,55,44,1)')  # Highlight the queried node
        else:
            node_colors.append('rgba(75,125,180,1)')

    # For edges - add deltamz_int attribute as hover info
    edge_hover_text = []
    for edge in cluster_edges:
        edge_attr = G.edges[edge]

        # Create edge hover text
        edge_text = f"<b>Edge: {edge[0]} → {edge[1]}</b><br>"
        edge_text += f"Δm/z: {edge_attr.get('deltamz_int', 'N/A')}<br>"

        # Add any other edge attributes you want to display
        if 'score' in edge_attr:
            edge_text += f"Cosine Score: {edge_attr['score']:.3f}<br>"
        if 'mass_diff' in edge_attr:
            edge_text += f"Mass Diff: {edge_attr['mass_diff']:.4f}<br>"

        edge_hover_text.append(edge_text)

    # Create edge trace
    edge_trace = go.Scatter(
        x=x_edges, y=y_edges,
        line=dict(width=2, color='gray'),
        hoverinfo='text',
        hovertext=edge_hover_text,
        mode='lines+text',
        showlegend=False
    )
    
    # Create node trace
    node_trace = go.Scatter(
        x=x_nodes, y=y_nodes,
        mode='markers+text',
        hoverinfo='text',
        hovertext=hover_text,
        marker=dict(
            size=30,
            color=node_colors,
            line=dict(width=2, color='white')
        ),
        text=[str(node) for node in cluster_nodes],
        textposition="middle center",
        textfont=dict(size=10, color='lightgrey'),
        showlegend=False
    )
    
    # Create figure
    fig = go.Figure(data=[edge_trace, node_trace])

    fig.update_layout(
        title=f'Component {component_id} - Contains Node {node_id}',
        showlegend=False,
        hovermode='x',
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='white',
        width=width,
        height=height
    )

    print(edge_hover_text)
    print(x_edges)
    print(y_edges)
    # Add cluster info
    identified = sum(1 for node in cluster_nodes
                    if 'library_compound_name' in G.nodes[node] 
                    and G.nodes[node]['library_compound_name'])
    
    info_text = f"Component {component_id}<br>"
    info_text += f"Nodes: {len(cluster_nodes)}<br>"
    info_text += f"Edges: {subgraph.number_of_edges()}<br>"
    info_text += f"Identified: {identified}<br>"
    info_text += f"<b>Queried node: {node_id}</b>"
    
    fig.add_annotation(
        x=0.02, y=0.98,
        xref='paper', yref='paper',
        text=info_text,
        showarrow=False,
        align='left',
        bgcolor='rgba(255,255,255,0.9)',
        # bordercolor='red',
        borderwidth=2,
        font=dict(size=16)
    )
    
    return fig

def get_cluster_id(G, node_id):
    """
    Simple function to just get the component ID for a node.
    
    Parameters:
    G: NetworkX graph
    node_id: The node ID to look up (same as Feature ID
    
    Returns:
    component_id: The component ID, or None if not found
    """
    if node_id not in G.nodes():
        # print(f"❌ Node {node_id} not found!")
        return None

    component_id = G.nodes[node_id].get('component')
    
    return component_id

def quick_cluster_plot(G, node_id):
    """
    One-liner function: provide node ID, get the plot.
    """
    fig = plot_cluster_by_node(G, node_id)
    if fig:
        fig.show()
    return fig

# # Example usage:
# """
# # Load your network
# G = nx.read_graphml('your_network.graphml')
#
# # Method 1: Just get the component ID
# component_id = get_component_id(G, 'your_node_id')
#
# # Method 2: Get component ID and plot the cluster
# fig = plot_cluster_by_node(G, 'your_node_id')
# fig.show()
#
# # Method 3: One-liner
# quick_cluster_plot(G, 'your_node_id')
# """
#
# print("🎯 Simple cluster plotter ready!")
# print("\nUsage:")
# print("1. component_id = get_component_id(G, 'node_123')")
# print("2. fig = plot_cluster_by_node(G, 'node_123'); fig.show()")
# print("3. quick_cluster_plot(G, 'node_123')  # One-liner")