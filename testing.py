import numpy as np
from Algebraic_equivalence import model_contained
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def compare_graphs(graph_list):
    n_graphs = len(graph_list)
    grid = np.empty((n_graphs, n_graphs), dtype="object")
    grid_numeric = np.empty((n_graphs, n_graphs))
    mapping = {True: 1, False: 0, "INC": 2}
    
    for graph1_index in tqdm(range(n_graphs)):
        for graph2_index in range (n_graphs):
            graph1 = graph_list[graph1_index]
            graph2 = graph_list[graph2_index]
            
            contained = model_contained(graph1, graph2, method="EID_TSID")
            
            grid[graph1_index, graph2_index] = contained
            grid_numeric[graph1_index, graph2_index] = mapping[contained]
            
    return grid_numeric


def plot_grid(grid_numeric, graph_ids=[], save=False):
    colors = ["#e82c39", "#4257f5", "#0a0000"]
    cmap = mcolors.ListedColormap(colors)

    fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    ax.imshow(grid_numeric, cmap=cmap, vmin=0, vmax=2)

    ax.set_xticks(np.arange(-0.5, np.shape(grid_numeric)[0], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, np.shape(grid_numeric)[0], 1), minor=True)
    ax.grid(which="minor", color="black", linestyle='-', linewidth=0.5)

    if graph_ids:  
        num_graphs = len(graph_ids)
        font_size = 10 - num_graphs/19
        ax.set_xticks(np.arange(num_graphs))
        ax.set_xticklabels(graph_ids, rotation=90, fontsize=font_size)
        
        ax.set_yticks(np.arange(len(graph_ids)))
        ax.set_yticklabels(graph_ids, fontsize=font_size)
    
    ax.set_xlabel("G'")
    ax.set_ylabel("G")

    if save == True:
        plt.savefig('plot.png', dpi=300)
        
    plt.show()

