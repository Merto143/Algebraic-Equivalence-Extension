import numpy as np
from Algebraic_equivalence import model_contained
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from DMG import DMG

def compare_graphs(graph_list: list["DMG"]):
    """
    Compare all pairs of graphs in the input list using the EID+TSID containment method.

    For each pair (G1, G2), check whether G1 is algebraically contained in G2. The result is stored
    in a numeric grid using the following mapping:
        - True  → 1  (G1 is contained in G2)
        - False → 0  (G1 is not contained in G2)
        - "INC" → 2  (Containment check was inconclusive)

    Args:
        graph_list (list[DMG]): A list of directed mixed graphs to compare.

    Returns:
        np.ndarray: A numeric matrix where entry (i, j) corresponds to the result of
                    model_contained(graph_list[i], graph_list[j]).
    """
    n_graphs = len(graph_list)

    # Initialize result grid
    grid_numeric = np.empty((n_graphs, n_graphs))                       # Numeric results (1, 0, 2)
    mapping = {True: 1, False: 0, "INC": 2}                             # Define mapping from symbolic to numeric results
    
    # Compare each pair of graphs
    for graph1_index in tqdm(range(n_graphs)):
        for graph2_index in range (n_graphs):
            graph1 = graph_list[graph1_index]
            graph2 = graph_list[graph2_index]
            
            # Check if graph1 is contained in graph2 using the EID+TSID method
            contained = model_contained(graph1, graph2, method="EID_TSID")
            
            # Store symbolic and numeric results
            grid_numeric[graph1_index, graph2_index] = mapping[contained]
            
    return grid_numeric


def plot_grid(grid_numeric:  np.ndarray[int], graph_ids: list[tuple[int, int]] = [], save: bool = False):
    """
    Plot a heatmap of the numeric model containment grid.

    Each cell (i, j) of the grid indicates whether graph Gi is contained in Gj using:
        - 0 (red): Not contained
        - 1 (blue): Contained
        - 2 (black): Inconclusive

    Args:
        grid_numeric (np.ndarray): A 2D numeric grid with containment results (0, 1, or 2).
        graph_ids (list, optional): A list of graph names or labels to use as axis ticks.
        save (bool, optional): If True, saves the plot as 'plot.png'.
    """
    # Define a custom colormap: 0 → red, 1 → blue, 2 → black
    colors = ["#e82c39", "#4257f5", "#0a0000"]
    cmap = mcolors.ListedColormap(colors)

    # Create a square figure and plot the grid as an image
    fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    ax.imshow(grid_numeric, cmap=cmap, vmin=0, vmax=2)

     # Add minor ticks between grid cells to draw dividing lines
    ax.set_xticks(np.arange(-0.5, np.shape(grid_numeric)[0], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, np.shape(grid_numeric)[0], 1), minor=True)
    ax.grid(which="minor", color="black", linestyle='-', linewidth=0.5)

    # Set custom axis labels if graph IDs are provided
    if graph_ids:  
        num_graphs = len(graph_ids)
        font_size = 10 - num_graphs/19
        ax.set_xticks(np.arange(num_graphs))
        ax.set_xticklabels(graph_ids, rotation=90, fontsize=font_size)
        
        ax.set_yticks(np.arange(len(graph_ids)))
        ax.set_yticklabels(graph_ids, fontsize=font_size)
    
    ax.set_xlabel("G'")
    ax.set_ylabel("G")

     # Save the plot to file if requested
    if save == True:
        plt.savefig('plot.png', dpi=300)
        
    # Show the plot
    plt.show()

