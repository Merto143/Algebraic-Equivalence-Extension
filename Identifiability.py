from treks import htr, tr
from flow import max_flow, max_flow_TSID
import numpy as np
from itertools import chain, combinations


def HTC_identifiable(graph):
    graph.reset_id_grid()
    HTC_order = []
    var_order = []
    change = True
    
    for var in graph.variables:
        if len(var.parents) == 0:
            HTC_order.append([var, []])  
            var_order.append(var)
                 

    while change == True and len(HTC_order) != len(graph.variables):
        change = False
        
        for v in graph.variables:
            if not v in var_order:
                htr_v = htr(graph, v)
                v_min_htr = set(graph.variables) - set(htr_v)
                A = (v_min_htr | set(var_order)) - set([v]) - set(v.siblings)

                flow, flow_dict = max_flow(graph, v, A)          

                if flow == len(v.parents):
                    # Save the set of HTC identifiable nodes 
                    Yv = []
                    for node in flow_dict["s"]:
                        if flow_dict["s"][node] == 1:
                            Yv.append(graph.str_to_var[node[2]])

                    # Add variable to solved nodes
                    var_order.append(v)
                    HTC_order.append([v, Yv])
                    v_id = graph.var_to_int[v]

                    for parent in v.parents:
                        parent_id = graph.var_to_int[parent]
                        graph.id_grid[parent_id, v_id] = 1

                    change = True

    HTC_identifiable = len(HTC_order) == len(graph.variables)
    
    if HTC_identifiable:
        return True, HTC_order
    
    else:
        return False, None
    

def EID_step(graph):
    change = False
    order = []

    for var_id, var in enumerate(graph.variables):
        unsolved = []
        maybe = []

        # obtain UNSOLVED set
        for w_id, w in enumerate(graph.variables):
            if graph.id_grid[w_id, var_id] == -1:
                unsolved.append(w)
        
        # obtain MAYBE_ALLOWED set
        Y = set(graph.variables) - {var} - set(var.siblings)
        
        for y in Y: 
            Z = set(htr(graph, var)) & set(y.parents)
            known = True
                
            for z in Z:
                
                if graph.id_grid[graph.var_to_int[z], graph.var_to_int[y]] == -1:
                    known = False

            if known == True:
                maybe.append(y)
        
        combs = list(chain.from_iterable(combinations(unsolved, r) for r in range(1, len(unsolved) + 1)))

        allowed = []
        
        for W in combs:
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

            flow, flow_dict = max_flow(graph, var, allowed)
            
            if len(W) == flow:
                change = True
                Y = []

                for node in flow_dict["s"]:
                        if flow_dict["s"][node] == 1:
                            Y.append(graph.str_to_var[node[2]])
                
                for w in W:                
                    w_id = graph.var_to_int[w]  
                    graph.id_grid[w_id, var_id] = 1

                Y = sorted(Y)
                W = sorted(list(W))

                order.append([var, W, Y])


            if change == True:
                break
                            
    return change, order


def EID_identifiable(graph):
    graph.reset_id_grid()
    change = True

    EID_order = []

    while change == True:
        change, order = EID_step(graph)

        if change:
            EID_order = EID_order + order

    if not np.any(graph.id_grid == -1):
        identifiable = True
    else:
        identifiable = False

    return identifiable, EID_order


def TSID_step(graph):
    change = False

    for w0_id, w0 in enumerate(graph.variables):
        for var_id, var in enumerate(graph.variables):
            if graph.id_grid[w0_id, var_id] == -1:
                remove_edges = [[w0, var]]
                
                for parent in var.parents:
                    parent_id = graph.var_to_int[parent]
                    
                    if graph.id_grid[parent_id, var_id] == 1:
                        remove_edges.append([parent, var])

                T_set = set(graph.variables) - {w0} - {var}
                T_combs = list(chain.from_iterable(combinations(T_set, r) for r in range(1, graph.N_nodes - 1)))
                
                for T_comb in T_combs:
                    if not set(var.descendents) & (set(T_comb) | {var}):
                        S_combs = list(combinations(graph.variables, len(T_comb) + 1))

                        for S_comb in S_combs:
                            max_flow, flow_dict = max_flow_TSID(graph, S_comb, set(T_comb) | {w0})

                            if max_flow == len(S_comb):
                                max_flow_prime, flow_dict_prime = max_flow_TSID(graph, S_comb, set(T_comb) | {var}, remove_edges=remove_edges)
                                if  max_flow_prime < len(S_comb):
            
                                    graph.id_grid[w0_id, var_id] = 1
                                    change = True

                                    S_comb = sorted(S_comb)
                                    T_comb = sorted(T_comb)
                                                             
                                    return change, [w0, var, S_comb, T_comb]
                                    
                                                        
    return change, []

def TSID_identifiable(graph):
    graph.reset_id_grid()
    change = True
    TSID_order = []

    while change:
        change, order = TSID_step(graph)

        if change:
            TSID_order.append(order)

    if not np.any(graph.id_grid == -1):
        identifiable = True
    else:
        identifiable = False

    return identifiable, TSID_order


def EID_TSID_identifiable(graph): 
    graph.reset_id_grid()
    total_change = True
    total_order = []

    while total_change:
        total_change = False
        EID_change = True
        TSID_change = True
        EID_order = []
        TSID_order = []

        while EID_change:
            EID_change, order = EID_step(graph)

            if EID_change:
                EID_order = EID_order + order
                total_change = True
    
        while TSID_change:
            TSID_change, order = TSID_step(graph)

            if TSID_change:
                TSID_order.append(order)
                total_change = True
        
        if EID_order:
            total_order.append(["EID", EID_order])  

        if TSID_order:
            total_order.append(["TSID", TSID_order])

    if not np.any(graph.id_grid == -1):
        identifiable = True
    else:
        identifiable = False

    return identifiable, total_order

        
        