from treks import htr, tr
from flow import max_flow, max_flow_TSID
import numpy as np
from itertools import chain, combinations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from DMG import DMG 

def HTC_identifiable(graph: "DMG"):
    """
        Use the half-trek criterion Foygel et al. (2012) to check whether a DMG is HTC-identifiable 

        Args:
            graph (DMG): The DMG to analyze.

        Returns:
            HTC_identifiable (bool): True if the DMG is HTC-identifiable, False otherwise
            HTC_order (list[Tuple[Variable, list[Variable]]]): The HTC order if the graph is HTC identifiable.
            Each entry is a pair (v, Yv), where:
                - v is a node in the graph.
                - Yv is the HTC-identifying set for v.
    """
    graph.reset_id_grid()                           # Reset id_grid (used to track which edges are identified)
    HTC_order = []
    var_order = []
    change = True                                   # Track if new variables were identified in the current pass
    
    # Identify all nodes with no parents
    for var in graph.variables:
        if len(var.parents) == 0:
            HTC_order.append([var, []])             # No identifying set needed
            var_order.append(var)                   # Mark as identified
                 

    while change and len(HTC_order) != len(graph.variables):
        change = False                              
        
        # Go through all non identified variables. 
        for v in graph.variables:
            if not v in var_order:
                htr_v = htr(graph, v)               # Obtain all half-trek reachable nodes from v
                v_min_htr = set(graph.variables) - set(htr_v)   
                A = (v_min_htr | set(var_order)) - set([v]) - set(v.siblings)  # Obtain the set of allowed nodes

                flow, flow_dict = max_flow(graph, v, A)  # Construct the flow-graph and calculate maximum flow
                
                # The node is identifiable if the max flow is equal to the number of parents of v
                if flow == len(v.parents):
                    # Save the set of HTC identifiable nodes 
                    Yv = []  
                    for node in flow_dict["s"]:
                        if flow_dict["s"][node] == 1:
                            Yv.append(graph.str_to_var[node[2]])

                    var_order.append(v)  # Mark node as identified
                    HTC_order.append([v, Yv])  
                    
                    # Mark edges from pa(v) to v as identified
                    v_id = graph.var_to_int[v]
                    for parent in v.parents:
                        parent_id = graph.var_to_int[parent]
                        graph.id_grid[parent_id, v_id] = 1   # Mark edge as identified in the ID grid

                    change = True

    # Check whether all nodes have been identified. 
    HTC_identifiable = len(HTC_order) == len(graph.variables)
    
    if HTC_identifiable:
        return True, HTC_order
    
    else:
        return False, None
    

def EID_step(graph: "DMG"):
    """
        Perform one iteration of the edgewise identificaton algorithm on a DMG.

        This function attempts to identify edges in the DMG using the EID criterion as described by Weihs et al (2018).

        Args:
            graph (DMG): The DMG to analyze.

        Returns:
            change (bool): True if the new edges have been identified in this sweep, False otherwise
            order (list[Tuple[Variable, list[Variable]]]): A list of identification steps performed. 
            Each entry in order is a triple (v, W, Y), where:
                - v is a node in the graph.
                - W is a subset of parents nodes of v
                - Y is the EID-identifying set corresponding to W for v.
    """
    change = False                                                  # Track if new edges are identified in the current step
    order = []

    # Go once through every variable
    for var_id, var in enumerate(graph.variables):
        unsolved = []
        maybe = []

        # Build UNSOLVED set
        for w_id, w in enumerate(graph.variables):
            if graph.id_grid[w_id, var_id] == -1:
                unsolved.append(w)                                  # List of parents of v that are not yet identified
        
        # Construct MAYBE_ALLOWED set
        Y = set(graph.variables) - {var} - set(var.siblings)        
        for y in Y: 
            Z = set(htr(graph, var)) & set(y.parents)               # Calculate htr(v) intersection tr(pa(v))
            known = True
            
            # Check if all z -> y are identified
            for z in Z:
                if graph.id_grid[graph.var_to_int[z], graph.var_to_int[y]] == -1:
                    known = False

            if known == True:
                maybe.append(y)
        
        # Check all possible subsets of  UNSOLVED as potential W sets
        combs = list(chain.from_iterable(combinations(unsolved, r) for r in range(1, len(unsolved) + 1)))

        allowed = []
        
        # Try each subset W of unsolved parents
        for W in combs:
            # For each W, construct a set of allowed nodes
            for y in maybe:
                tr_set = []
                
                for p in y.parents:
                    if graph.id_grid[graph.var_to_int[p], graph.var_to_int[y]] == -1:
                        tr_set = tr_set + tr(graph, p)

                tr_set = set(tr_set)
                htr_set = set(htr(graph, y))
                combined_set = (tr_set | htr_set) & set(unsolved)
                
                if combined_set.issubset(set(W)) and combined_set:
                    allowed.append(y)

            flow, flow_dict = max_flow(graph, var, allowed)                                 # Calculate the maxflow from allowed to var
            
            # If maxflow is equal to |W|, the identification is succesful
            if len(W) == flow:
                change = True
                Y = []

                # Obtain EID identifying set Y from flow output
                for node in flow_dict["s"]:
                        if flow_dict["s"][node] == 1:
                            Y.append(graph.str_to_var[node[2]])
                
                # Mark edges w -> var as identified
                for w in W:                
                    w_id = graph.var_to_int[w]  
                    graph.id_grid[w_id, var_id] = 1                                         # Mark edge as identified in the ID grid

                # Sort for consistent output order
                Y = sorted(Y)
                W = sorted(list(W))

                # Add EID step to the identification order
                order.append([var, W, Y])

            if change:
                break
                            
    return change, order


def EID_identifiable(graph: "DMG"):
    """
        Implementation of the edgewise identification algorithm (Weihs et al. 2018).
        This function checks whether a DMG is edgewise identifiable. 

        Args:
            graph (DMG): The DMG to analyze.

        Returns:
            EID_identifiable (bool): True if the DMG is EID-identifiable, False otherwise
            EID_order (list[Tuple[Variable, list[Variable], list[Variable]]]): The EID order if the graph is EID identifiable.
            Each entry is a triple (v, W, Y), where:
                - v is a node in the graph.
                - W is a subset of parents nodes of v
                - Y is the EID-identifying set corresponding to W for v.
    """
    graph.reset_id_grid()                        # Reset id_grid (used to track which edges are identified)
    change = True                                # Track if new edges were identified in the current pass
    EID_order = []                               

    # Apply EID steps untill no more edges can be identified
    while change:
        change, order = EID_step(graph)         # Attempt to identify more edges

        if change:
            EID_order = EID_order + order       # Add newly identified entry to the full order

    # Check whether all edges are identified using the ID grid
    if not np.any(graph.id_grid == -1):
        EID_identifiable = True
    else:
        EID_identifiable = False

    return EID_identifiable, EID_order


def TSID_step(graph: "DMG"):
    """
        Perform one iteration of the trek-separation identificaton algorithm on a DMG.

        This function attempts to identify individuel using the TSID criterion as described by Weihs et al (2018).

        Args:
            graph (DMG): The DMG to analyze.

        Returns:
            change (bool): True if a new edge has been identified in the current step, False otherwise.
            order: The TSID identifying variables of the edge that has been identified.
            The entry is a quadruple (w0, v, S, T), where:
                - v is a node in the graph.
                - w0 is a parent of v. 
                - S is a subset of source nodes.
                - T is a subset of sink nodes.
    """
    change = False                                              # Track if new edges were identified in the current pass

     # Try identifying each edge w0 → v that is currently unidentifiable
    for w0_id, w0 in enumerate(graph.variables):
        for var_id, var in enumerate(graph.variables):
            if graph.id_grid[w0_id, var_id] == -1:
                remove_edges = [[w0, var]]                      # Edge to test for identifiability
                
                 # Also remove other already-identified incoming edges to v
                for parent in var.parents:
                    parent_id = graph.var_to_int[parent]
                    
                    if graph.id_grid[parent_id, var_id] == 1:
                        remove_edges.append([parent, var])

                # Build all candidate sink sets T 
                T_set = set(graph.variables) - {w0} - {var}
                T_combs = list(chain.from_iterable(combinations(T_set, r) for r in range(1, graph.N_nodes - 1)))
                
                for T_comb in T_combs:
                    # Skip T if v has descendants in T ∪ {v} — violates TSID condition
                    if not set(var.descendants) & (set(T_comb) | {var}):
                        S_combs = list(combinations(graph.variables, len(T_comb) + 1))              # Source set S must be larger than T — generate S accordingly

                        for S_comb in S_combs:
                            # Compute flow from S to T ∪ {w0} (edge present)
                            max_flow, flow_dict = max_flow_TSID(graph, S_comb, set(T_comb) | {w0})

                            if max_flow == len(S_comb):
                                max_flow_prime, flow_dict_prime = max_flow_TSID(graph, S_comb, set(T_comb) | {var}, remove_edges=remove_edges)

                                # If flow is strictly smaller, w0 → v is identifiable
                                if  max_flow_prime < len(S_comb):
            
                                    graph.id_grid[w0_id, var_id] = 1                                # Mark edge as identified
                                    change = True

                                    # Sort for consistant output.
                                    S_comb = sorted(S_comb)
                                    T_comb = sorted(T_comb)
                                    order = [w0, var, S_comb, T_comb]
                                                             
                                    return change, order
                                    
    order = []  # No edges identified in this step
    return change, order

def TSID_identifiable(graph: "DMG"):
    """
        Implementation of the trek-saperation identification algorithm (Weihs et al. 2018).
        This function checks whether a DMG is trek-separation identifiable. 

        Args:
            graph (DMG): The DMG to analyze.

        Returns:
            TSID_identifiable (bool): True if the DMG is TSID-identifiable, False otherwise
            TSID_order (list[Tuple[Variable, list[Variable], list[Variable]]]): The TSID order if the graph is TSID identifiable.
            Each entry is a quadruple (w0, v S, T), where:
                - v is a node in the graph.
                - w0 is a parent of v. 
                - S is a subset of source nodes.
                - T is a subset of sink nodes.
    """
    graph.reset_id_grid()                           # Reset id_grid (used to track which edges are identified)
    change = True                                   # Track if new edges were identified in the current pass
    TSID_order = []

    # Apply TSID steps untill no more edges can be identified
    while change:
        change, order = TSID_step(graph)            # Attempt to identify more edges

        if change:
            TSID_order.append(order)                # Add newly identified entry to the full order

    # Check whether all edges are identified using the ID grid
    if not np.any(graph.id_grid == -1):
        TSID_identifiable = True
    else:
        TSID_identifiable = False

    return TSID_identifiable, TSID_order


def EID_TSID_identifiable(graph: "DMG"): 
    """
        Implementation of the EID+TSID identification algorithm.

        This function checks whether a DMG is identifiable using a combination of edgewise identification and trek-separation identification,
        as defined by (Weihs et al. 2018).

        Args:
            graph (DMG): The DMG to analyze.

        Returns:
            EID_TSID_identifiable (bool): True if the DMG is EID+TSID-identifiable, False otherwise
            total_order: The identification order if the graph is EID+TSID-identifiable.
            Each step in total_order is a tuple (method, order), where:
                - method (str): "EID" edgewise identification was used in that step
                                or "TSID" trek-separation identification was used
                - order: The EID or TSID order used in that step.
    """
    graph.reset_id_grid()                                       # Reset id_grid (used to track which edges are identified)
    total_change = True                                         # Track if new edges were identified in the current pass
    total_order = []

    # Iteratlively apply EID and TSID untill no more edges can be identified
    while total_change:
        total_change = False
        EID_change = True
        TSID_change = True
        EID_order = []                                          # Save identification order of steps that use EID
        TSID_order = []                                         # Save identification order of steps that use TSID

        # Apply EID untill no more new edges can be identified
        while EID_change:
            EID_change, order = EID_step(graph) 

            if EID_change:
                EID_order = EID_order + order
                total_change = True
    
        # Apply TSID untill no more new edges can be identified
        while TSID_change:
            TSID_change, order = TSID_step(graph)

            if TSID_change:
                TSID_order.append(order)
                total_change = True
        
        # Append EID order to the total order
        if EID_order:
            total_order.append(["EID", EID_order])  

        # Append TSID order to the total order
        if TSID_order:
            total_order.append(["TSID", TSID_order])

    # Check whether all edges are identified using the ID grid
    if not np.any(graph.id_grid == -1): 
        EID_TSID_identifiable = True
    else:
        EID_TSID_identifiable = False

    return EID_TSID_identifiable, total_order

        
        