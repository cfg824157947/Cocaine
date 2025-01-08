import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def visualize_self_loops():
    # Node position
    x, y = 0.5, 0.5
    node_radius = 0.1

    # Loop parameters
    loop_radius = node_radius/2# Radius of the self-loop
    loop_offset = 0.05  # Offset from the node's center

    # Parametric equation for the loop
    theta = np.linspace(0, 2 * np.pi, 100)
    loop_center = (x, y + loop_offset)  # Offset the loop vertically
    loop_x = loop_center[0] + loop_radius * np.cos(theta)
    loop_y = loop_center[1] + loop_radius * np.sin(theta)

    # Plot the node
    plt.scatter(x, y, s=2000, color='skyblue', zorder=2)
    plt.text(x, y, 'NodeNode', fontsize=12, ha='center', va='center', zorder=3)

    # Plot the self-loop
    plt.plot(loop_x, loop_y, color='red', zorder=1)

    # Add arrowhead for the self-loop
    arrow_pos = loop_center[0] + loop_radius * np.cos(np.pi / 4), loop_center[1] + loop_radius * np.sin(np.pi / 4)
    plt.arrow(
        arrow_pos[0], arrow_pos[1], 0.01, 0.01,  # Small arrow to indicate direction
        head_width=0.03, head_length=0.05, fc='red', ec='red'
    )

    # Adjust plot
    plt.gca().set_aspect('equal')
    plt.axis('off')
    plt.title("Visualization of a Self-Loop Arrow")
    save_path = "self_loop.png"
    plt.savefig(save_path, dpi=300, bbox_inches='tight', format=save_path.split('.')[-1])
    plt.show()

plt.clf()
visualize_self_loops()




def plot_network_with_save(edges_df, node_size=2000, save_path=None, dpi=300):
    """
    Plot a directed network graph with customizable node sizes and save the figure.

    Parameters:
        edges_df (DataFrame): DataFrame containing edges with two columns (source, target).
        node_size (float): Size of the nodes in scatter (affects visual size).
        save_path (str): Path to save the figure (e.g., 'figure.png'). If None, the figure is not saved.
        dpi (int): Dots per inch for the saved figure.
    """
    # Assign positions to nodes in a circular layout
    nodes = pd.concat([edges_df[0], edges_df[1]]).unique()
    num_nodes = len(nodes)
    angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint=False)
    positions = {node: (np.cos(angle), np.sin(angle)) for node, angle in zip(nodes, angles)}
    node_radius = 0.1  # Adjust radius as needed

    # Create the figure
    plt.figure(figsize=(10, 10))

    # Plot arrows for edges
    for _, row in edges_df.iterrows():
        start, end = row[0], row[1]
        start_pos = np.array(positions[start])
        end_pos = np.array(positions[end])

        # Calculate direction and adjust start and end positions
        direction = end_pos - start_pos
        length = np.linalg.norm(direction)
        if length != 0:  # Avoid division by zero for self-loops
            unit_direction = direction / length
            start_adjusted = start_pos + node_radius * unit_direction
            end_adjusted = end_pos - node_radius * unit_direction
            #start_adjusted = start_pos + node_radius * unit_direction
            #end_adjusted = end_pos - node_radius * unit_direction

            # Plot the arrow
            plt.arrow(
                start_adjusted[0], start_adjusted[1],
                *(end_adjusted - start_adjusted),
                head_width=0.05, head_length=0.1,
                fc='gray', ec='gray', linewidth=1.0,
                length_includes_head=True
            )

    # Plot nodes with customizable size
    for node, (x, y) in positions.items():
        plt.scatter(x, y, s=node_size, color='skyblue', zorder=2)
        plt.text(x, y, str(node), fontsize=12, ha='center', va='center', zorder=3)

    # Adjust plot
    plt.gca().set_aspect('equal')
#    plt.axis('off')
    plt.title("Directed Network Graph")

    # Save the figure if save_path is specified
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', format=save_path.split('.')[-1])
        print(f"Figure saved to {save_path}")

    # Show the plot
    plt.show()
plot_network_with_save(edges_df, save_path=save_path)
plt.clf()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_network_with_groups(edges_df, node_size=300, save_path=None, dpi=300):
    """
    Plot a directed network graph with groups and distinct colors for nodes.
    
    Parameters:
        edges_df (DataFrame): DataFrame containing edges with two columns (source, target).
        node_size (float): Size of the nodes in scatter (affects visual size).
        save_path (str): Path to save the figure (e.g., 'figure.png'). If None, the figure is not saved.
        dpi (int): Dots per inch for the saved figure.
    """
    # Define node groups
    nodes = pd.concat([edges_df[0], edges_df[1]]).unique()
    group_colors = {
        "yellow": (0, 25),
        "green": (25, 50),
        "blue": (50, 75),
        "red": (75, 100)
    }
    node_colors = {}
    positions = {}

    # Assign positions to nodes in their respective groups
    plt.figure(figsize=(8, 8))
    num_groups = len(group_colors)
    group_angle_span = 2 * np.pi / num_groups  # Divide space equally

    for i, (color, (start, end)) in enumerate(group_colors.items()):
        group_nodes = [node for node in nodes if start <= int(node) < end]
        num_nodes = len(group_nodes)
        angles = np.linspace(
            i * group_angle_span,
            (i + 1) * group_angle_span,
            num_nodes,
            endpoint=False
        )
        for node, angle in zip(group_nodes, angles):
            positions[node] = (np.cos(angle), np.sin(angle))
            node_colors[node] = color

    # Plot edges
    for _, row in edges_df.iterrows():
        start, end = row[0], row[1]
        start_pos = np.array(positions[start])
        end_pos = np.array(positions[end])

        # Calculate direction and adjust start and end positions
        direction = end_pos - start_pos
        length = np.linalg.norm(direction)
        if length != 0:  # Avoid division by zero for self-loops
            unit_direction = direction / length
            node_radius = 0.15
            start_adjusted = start_pos + node_radius * unit_direction
            end_adjusted = end_pos - node_radius * unit_direction

            # Plot the arrow
            plt.arrow(
                start_adjusted[0], start_adjusted[1],
                *(end_adjusted - start_adjusted),
                head_width=0.05, head_length=0.1,
                fc='gray', ec='gray', linewidth=1.0,
                length_includes_head=True
            )

    # Plot nodes with their colors
    for node, (x, y) in positions.items():
        plt.scatter(x, y, s=node_size, color=node_colors[node], zorder=2)
        plt.text(x, y, str(node), fontsize=10, ha='center', va='center', zorder=3)

    # Adjust plot
    plt.gca().set_aspect('equal')
    plt.axis('off')
    plt.title("Directed Network Graph with Node Groups")

    # Save the figure if save_path is specified
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', format=save_path.split('.')[-1])
        print(f"Figure saved to {save_path}")

    # Show the plot
    plt.show()


# Example usage
# Generate edges DataFrame
edges_data = {
    0: np.random.randint(0, 100, 200),  # Source nodes
    1: np.random.randint(0, 100, 200)   # Target nodes
}
edges_df = pd.DataFrame(edges_data)

# Plot the network
plot_network_with_groups(edges_df, save_path="grouped_network_graph.png")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_network_with_rectangular_groups(edges_df, node_size=2000, save_path=None, dpi=300):
    """
    Plot a directed network graph with groups and distinct colors for nodes in a rectangular layout.

    Parameters:
        edges_df (DataFrame): DataFrame containing edges with two columns (source, target).
        node_size (float): Size of the nodes in scatter (affects visual size).
        save_path (str): Path to save the figure (e.g., 'figure.png'). If None, the figure is not saved.
        dpi (int): Dots per inch for the saved figure.
    """
    # Define node groups
    nodes = pd.concat([edges_df[0], edges_df[1]]).unique()
    group_colors = {
        "yellow": (0, 25),
        "green": (25, 50),
        "blue": (50, 75),
        "red": (75, 100)
    }
    node_colors = {}
    positions = {}

    # Layout parameters
    group_regions = {
        "yellow": (-1, 1, 0, 1),  # Top-left quadrant
        "green": (0, 1, 0, 1),   # Top-right quadrant
        "blue": (-1, 1, -1, 0),  # Bottom-left quadrant
        "red": (0, 1, -1, 0)     # Bottom-right quadrant
    }

    # Layout parameters
    group_regions = {
        "yellow": (-1, 0, 0, 1),  # Top-left quadrant
        "green": (0, 1, 0, 1),   # Top-right quadrant
        "blue": (-1, 0, -1, 0),  # Bottom-left quadrant
        "red": (0, 1, -1, 0)     # Bottom-right quadrant
    }
#np.random.seed(2)
    # Assign positions to nodes in their respective groups
    for color, (start, end) in group_colors.items():
        group_nodes = [node for node in nodes if start <= int(node) < end]
        x_min, x_max, y_min, y_max = group_regions[color]
        
        # Define the center of the group region
        center_x = (x_min + x_max) / 2
        center_y = (y_min + y_max) / 2

        # Generate positions with a normal distribution
        std_dev = 0.4  # Adjust this value to control spread (smaller = tighter cluster)
        num_nodes = len(group_nodes)

        # Generate random positions within the specified region
        xs = np.random.uniform(x_min, x_max, num_nodes)
        ys = np.random.uniform(y_min, y_max, num_nodes)
        #xs = np.random.normal(loc=center_x, scale=std_dev, size=len(group_nodes))
        #ys = np.random.normal(loc=center_y, scale=std_dev, size=len(group_nodes))
        for node, x, y in zip(group_nodes, xs, ys):
            positions[node] = (x, y)
            node_colors[node] = color

    # Create the figure
    plt.figure(figsize=(15, 15))

    # Plot edges
    for _, row in edges_df.iterrows():
        start, end = row[0], row[1]
        start_pos = np.array(positions[start])
        end_pos = np.array(positions[end])

        # Calculate direction and adjust start and end positions
        direction = end_pos - start_pos
        length = np.linalg.norm(direction)
        if length != 0:  # Avoid division by zero for self-loops
            unit_direction = direction / length
            node_radius = 0.05
            start_adjusted = start_pos + node_radius * unit_direction
            end_adjusted = end_pos - node_radius * unit_direction

            # Plot the arrow
            plt.arrow(
                start_adjusted[0], start_adjusted[1],
                *(end_adjusted - start_adjusted),
                head_width=0.03, head_length=0.05,
                fc='gray', ec='gray', linewidth=1.0,
                length_includes_head=True
            )

    # Plot nodes with their colors
    for node, (x, y) in positions.items():
        plt.scatter(x, y, s=node_size, color=node_colors[node], zorder=2)
        plt.text(x, y, str(node), fontsize=10, ha='center', va='center', zorder=3)

    # Adjust plot
    plt.gca().set_aspect('equal')
    plt.axis('off')
    plt.title("Directed Network Graph with Rectangular Node Groups")

    # Save the figure if save_path is specified
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', format=save_path.split('.')[-1])
        print(f"Figure saved to {save_path}")

    # Show the plot
    plt.show()


# Example usage
# Generate edges DataFrame
edges_data = {
    0: np.random.randint(0, 100, 200),  # Source nodes
    1: np.random.randint(0, 100, 200)   # Target nodes
}
edges_df = pd.DataFrame(edges_data)

# Plot the network
plot_network_with_rectangular_groups(edges_df, save_path="rectangular_network_graph.png")

import plotly.graph_objects as go
import numpy as np
import pandas as pd

def plot_network_with_plotly(edges_df, node_size=20, save_path=None):
    """
    Plot a directed network graph with groups and distinct colors for nodes in a rectangular layout using Plotly.
    """
    # Define node groups
    nodes = pd.concat([edges_df[0], edges_df[1]]).unique()
    group_colors = {
        "yellow": (0, 25),
        "green": (25, 50),
        "blue": (50, 75),
        "red": (75, 100)
    }
    node_colors = {}
    positions = {}

    # Layout parameters
    group_regions = {
        "yellow": (-1, 0, 0, 1),  # Top-left quadrant
        "green": (0, 1, 0, 1),   # Top-right quadrant
        "blue": (-1, 0, -1, 0),  # Bottom-left quadrant
        "red": (0, 1, -1, 0)     # Bottom-right quadrant
    }

    # Assign positions to nodes in their respective groups
    for color, (start, end) in group_colors.items():
        group_nodes = [node for node in nodes if start <= int(node) < end]
        x_min, x_max, y_min, y_max = group_regions[color]

        # Grid dimensions
        grid_size = int(np.ceil(np.sqrt(len(group_nodes))))  # Square grid
        x_positions = np.linspace(x_min, x_max, grid_size)
        y_positions = np.linspace(y_min, y_max, grid_size)

        # Create a grid of positions
        grid_positions = [(x, y) for x in x_positions for y in y_positions]
        for node, pos in zip(group_nodes, grid_positions):
            positions[node] = pos
            node_colors[node] = color

    # Create Plotly figure
    fig = go.Figure()

    # Add edges as lines
    for _, row in edges_df.iterrows():
        start, end = row[0], row[1]
        start_pos = positions[start]
        end_pos = positions[end]

        fig.add_trace(go.Scatter(
            x=[start_pos[0], end_pos[0], None],
            y=[start_pos[1], end_pos[1], None],
            mode='lines',
            line=dict(color='gray', width=1),
            hoverinfo='none'
        ))

    # Add nodes as scatter points
    for node, (x, y) in positions.items():
        fig.add_trace(go.Scatter(
            x=[x],
            y=[y],
            mode='markers+text',
            marker=dict(
                size=node_size,
                color=node_colors[node],
                line=dict(width=1, color='black')
            ),
            text=str(node),
            textposition="top center",
            name=str(node)
        ))

    # Customize layout
    fig.update_layout(
        title="Directed Network Graph with Non-Overlapping Node Groups (Plotly)",
        showlegend=False,
        xaxis=dict(showgrid=False, zeroline=False),
        yaxis=dict(showgrid=False, zeroline=False),
        plot_bgcolor='white'
    )

    # Save to file
    if save_path:
        fig.write_image(save_path)
        print(f"Figure saved to {save_path}")

    # Show the plot
#    fig.show()


# Example usage
edges_data = {
    0: np.random.randint(0, 100, 200),  # Source nodes
    1: np.random.randint(0, 100, 200)   # Target nodes
}
edges_df = pd.DataFrame(edges_data)

plot_network_with_plotly(edges_df, save_path="plotly_network_graph.png")
