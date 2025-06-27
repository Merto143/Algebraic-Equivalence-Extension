from Variable import Variable

import networkx as nx
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from DMG import DMG
    from Variable import Variable


def max_flow(graph: "DMG", v: "Variable", A: list["Variable"]): 
    """
        Construct a flow graph for a DMG following Foygel et al. (2012) and compute the maximum flow

        Args:
            graph (DMG): The DMG for which to apply the max-flow algorithm.
            v (Variable): The target node.
            A: (List[Variable]): A conditioning set of variables.

        Returns:
            flow_value (int): The value of the maximum flow.
            flow_dict (dict): A dictionary describingt the flow on each edge.

    """
    g = nx.DiGraph()                                                        # Initialize a directed graph for the flow network
    g.add_edge("s", "t", capacity=0)                                        # Add a dummy source-sink edge

    # Step 1: 
    for var in A:            
        g.add_edge("s", f"l_{var.name}_in")                                 # Add s -> L(a) for each a in A
        g.add_edge(f"l_{var.name}_out", f"r_{var.name}_in")                 # Add L(a) -> R(a) for each a in A    
        for sibling in var.siblings:
            g.add_edge(f"l_{var.name}_out", f"r_{sibling.name}_in")         # Add L(a) -> R(w) for each a <-> w in B    

    # Step 2:
    for var in graph.variables:
        # Let all edges have max capacity 1
        g.add_edge(f"l_{var.name}_in", f"l_{var.name}_out", capacity=1)
        g.add_edge(f"r_{var.name}_in", f"r_{var.name}_out", capacity=1)
        
        for child in var.children: 
            g.add_edge(f"r_{var.name}_out", f"r_{child.name}_in")           # Add R(w) -> R(u) for each w -> u in D
    
    # Step 3:
    for parent in v.parents:
        g.add_edge(f"r_{parent.name}_out", "t")                             # Add R(w) -> t for each w in pa(v)

    flow_value, flow_dict = nx.maximum_flow(g, 's', 't')                    # Compute the maximum flow of the constructed flow graph
    
    return flow_value, flow_dict


def max_flow_TSID(graph: "DMG", S: list["Variable"], T: list["Variable"], remove_edges: list[list["Variable", "Variable"]] = []):
    """
        Construct a flow graph for a DMG following Weihs et al. (2018) and compute the maximum flow

        Args:
            graph (DMG): The DMG for which to apply the max-flow algorithm.
            S (List[Variable]): The source node set
            T (List[Variable]): The target node set
            removed_edges (List[List[Variable, Variable]]): A list of edges to remove when constructing G*

        Return:
            flow_value (int): The value of the maximum flow.
            flow_dict (dict): A dictionary describingt the flow on each edge.
    """
    g = nx.DiGraph()                                                        # Initialize a directed graph for the flow network     
    g.add_edge("s", "t", capacity=0)                                        # Add a dummy source-sink edge
    
    # Connect source nodes in S to the flow source 's'
    for var in S:            
        g.add_edge("s", f"l_{var.name}_in")

    # Connect target nodes in T to the flow sink 't'
    for var in T:
        g.add_edge(f"r_{var.name}_out", "t")
    
    for var in graph.variables:
        # Let all edges have max capacity 1
        g.add_edge(f"l_{var.name}_in", f"l_{var.name}_out", capacity=1)     
        g.add_edge(f"r_{var.name}_in", f"r_{var.name}_out", capacity=1)

        g.add_edge(f"l_{var.name}_out", f"r_{var.name}_in")                 # L(i) -> R(i) for each i in V

        for child in var.children:
            g.add_edge(f"l_{child.name}_out", f"l_{var.name}_in")           # L(i) -> L(j) for each j -> i in D
            g.add_edge(f"r_{var.name}_out", f"r_{child.name}_in")           # R(i) -> R(j) for each i -> j in D
        
        for sibling in var.siblings:
            g.add_edge(f"l_{var.name}_out", f"r_{sibling.name}_in")         # L(i) -> R(j) for each i <-> j in B
    
    # Remove edges form the flow graph as specified to create G*
    for parent, child in remove_edges:
        g.remove_edge(f"r_{parent.name}_out", f"r_{child.name}_in")

    flow_value, flow_dict = nx.maximum_flow(g, 's', 't')                    # Compute the maximum flow of the constructed flow graph

    
    return flow_value, flow_dict
    

