from DMG import DMG

from graphviz import Digraph
import itertools
import numpy as np
import string

from typing import Tuple, List
from numpy.typing import NDArray


def id_to_edges(nodes: List[str], graph_id: Tuple[int, int]):
    """
        Compute a list of directed and bidirected edges from a graph ID.

        Args:
            nodes (List[str]): A list of variable names.
            graph_id (Tuple[int, int]): a pair of bit-encoded integers, one for directed edges and one for bidirected edges.

        Return:
            edges (List[List[str, str]]): A list of directed and bidirected edges 
            omitted_nodes (List[str]): A list of omitted variable names.

    """
    omitted_nodes = []
    edges = []   
    omitted_id = 1                                      # Counter to assign unique names to omitted variables
    d, b = graph_id                                     # Save bit-encoded integers

    # Decode directed edges from bit-encoded integer
    # Iterate iver each ordered pair of nodes (excluding self-loops)
    for var1_id, var1 in enumerate(nodes):              
        for var2_id, var2 in enumerate(nodes):
            if not var1_id == var2_id:                  
                if d % 2 == 1:                          # If current bit is 1, add directed edge
                    edges.append((var1, var2))          
    
                d = d // 2                              # Shift to next bit

    # Decode bidirected edges from bit-encoded integer
    # Iterate over each unordered pair of nodes
    for var1_id, var1 in enumerate(nodes):              
        for var2_id, var2 in enumerate(nodes):
            if var2_id > var1_id:                       
                if b % 2 == 1:                          # If current bit is 1, add bidirected edge
                    o_var = f"o{omitted_id}"            # Uniquely name the new omitted variable
                    omitted_nodes.append(o_var)         # Add to list of omitted nodes
                    edges.append((o_var, var1))         # Add edge from omitted node to var1
                    edges.append((o_var, var2))         # Add edge from omitted node to var2
                    
                    omitted_id += 1                     
                    
                b = b // 2                              # Shift to the next bit
                
    return edges, omitted_nodes


def graph_to_id(graph: DMG):
    """
        Compute unique graph ID from a DMG.

        Args:
            graph (DMG): The graph from which to obtain the graph ID.

        Return:
            d (int): A bit encoded integer representing the directed edges.
            b (int): A bit encoded integer representing the bidirected edges.

    """
    variables = graph.variables[::-1]            # Reverse the variable list to math encoding order
    d = 0                                        # Encoding for directed edges
    b = 0                                        # Encoding for bidirected edges

    # Encode directed edges
    # Iterate over all ordered pairs (excluding self loops)
    for var1 in variables:      
        for  var2 in variables:
            if not var1 == var2:
                d *= 2                          # Shift bits to make room for the next edge

                if var2 in var1.children:       # If var1 -> var2 exists in the DMG
                    d += 1                      # Set current bit to 1

    # Encode bidirected edges
    # Iterate over all unordered pair of nodes
    for var1_id, var1 in enumerate(variables[1:]):
        for var2 in variables[:var1_id + 1]: 
            b *= 2                              # Shift bits to make room for the next edge
            if var2 in var1.siblings:           # If var1 <-> var2 exists in the DMG
                b += 1                          # Set current bit to 1
                
    return d, b


def permutate_graph(graph: DMG):
    """
        Create a list containing all possible permutations of a DMG.

        Args:
            graph (DMG): The input graph to permutate.

        Return:
            graph_list (List[DMG]): A list of DMGs, each corresponding to a different permutation.

    """
    graph.generate_parameters()                                                         # Generate a set of parameters for the original graph
    n = graph.N_nodes                                                                   # Number of nodes in the graph
    graph_list = []

    # Iterate over all permuations of node indeces
    for pi in itertools.permutations(list(range(n))):
        P = np.zeros(shape=(n,n), dtype=bool)                                           # Create permutation matrix P corresponding to pi
        for i in range(n):
            P[i, pi[i]] = True

        # Apply permutation to Lambda and Omega
        Lambda_new = P.dot(graph.Lambda).dot(P.T)
        Omega_new = P.dot(graph.Omega).dot(P.T)         

        graph_new = graph_from_parameters(Lambda_new, Omega_new, prime=graph.prime)     # Create graph from permutated parameters
        graph_list.append(graph_new)                                                    # Add permutated graph to the list

    return graph_list


def vizualize_graph(graph: DMG):
    """
        Create a vizual representation of a DMG using Graphviz.

        Args:
            graph (DMG): The input graph to vizualize.

        Return:
            g (Digraph): A Graphviz Digraph object representing the input graph.

    """
    g = Digraph('G', filename="graph.gv")                                           # Initialize Graphviz object
    added_pairs = []                                                                # Track which bidirected edges have been added to avoid duplicates
    
    for variable in graph.variables:
        g.node(variable.name)                                                       # Add a node for each variable  
        
        #Add directed edges
        for parent in variable.parents:     
            g.edge(parent.name, variable.name)                                      

        # Add bidirected edges
        for sibling in variable.siblings:
            if not (variable, sibling) in added_pairs:
                g.edge(variable.name, sibling.name, dir='both', style="dotted")
                added_pairs.append((sibling, variable))                             # Save bidirected edge to avoid duplicates

    return g


def graph_from_parameters(Lambda: NDArray, Omega: NDArray, prime: int = 2**31-1):
    """
        Construct a DMG from a given set of parameters.

        Args:
            Lambda (NDArray): A square matrix representing the directed edges.
            Omega (NDArray): A symmetric matrix representing the bidirected edges. 
            prime (int): The prime number defining the finite field for the DMG

        Return:
            graph (DMG): A DMG consistent with the structure of Lambda and Omega.
    """
    # Generate a list of unique variable names
    letters = string.ascii_lowercase
    num_variables = Lambda.shape[0]
    nodes = []

    for i in range(num_variables):
        name = ""
        n = i
        while True:
            n, rem = divmod(n, 26)
            name = letters[rem] + name
            if n == 0:
                break
            n -= 1                                                      # Adjust for 0-based indexing
        nodes.append(name)

    edges = []
    omitted_nodes = []
    omitted_count = 0                                                   # Counter to generate unique names for omitted variables
    
    # Add edges to the graph based on Lambda and Omega
    for row_id in range(num_variables):     
        for col_id in range(num_variables):
            # Add directed edges
            if Lambda[row_id, col_id]:
                edges.append((nodes[row_id], nodes[col_id]))        

            # Add bidirected edges
            if col_id > row_id:
                if Omega[row_id, col_id]:
                    # Generate unique omitted variable name
                    omitted_count += 1
                    omitted_nodes.append(f"o{omitted_count}")           
                    edges.append((f"o{omitted_count}", nodes[row_id]))
                    edges.append((f"o{omitted_count}", nodes[col_id]))

    graph = DMG(nodes, edges, omitted_nodes, prime=prime)               # Construct DMG from obtained nodes and edges

    return graph