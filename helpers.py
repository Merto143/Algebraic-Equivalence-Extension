from graphviz import Digraph
from DMG import DMG
import itertools
import numpy as np

def id_to_edges(nodes, graph_id):
    omitted_nodes = []
    edges = []   
    omitted_id = 1
    
    d_edges, b_edges = graph_id

    for var1_id, var1 in enumerate(nodes):
        for var2_id, var2 in enumerate(nodes):
            

            if not var1_id == var2_id:
                if d_edges % 2 == 1:
                    edges.append((var1, var2))
    
                d_edges = d_edges // 2
    
            if var2_id > var1_id:
                if b_edges % 2 == 1:
                    o_var = f"o{omitted_id}"
                    omitted_nodes.append(o_var)
                    edges.append((o_var, var1))
                    edges.append((o_var, var2))
                    
                    omitted_id += 1
                    
                b_edges = b_edges // 2
                
    return edges, omitted_nodes


def graph_to_id(graph):
    variables = graph.variables[::-1]
    d_id = 0
    b_id = 0
    for var1 in variables:
        for  var2 in variables:
            if not var1 == var2:
                d_id = d_id * 2

                if var2 in var1.children:
                    d_id += 1

    for var1_id, var1 in enumerate(variables[1:]):
        for var2 in variables[:var1_id+1]:
            b_id = b_id * 2
            if var2 in var1.siblings:
                b_id += 1
                
    return d_id, b_id


def permutate_graph(graph):
    graph.generate_parameters()
    d = graph.N_nodes

    graph_list = []
    for pi in itertools.permutations(list(range(d))):
        P = np.zeros(shape=(d,d), dtype=bool)
        for i in range(d):
            P[i, pi[i]] = True
            
        Lambda_new = P.dot(graph.Lambda).dot(P.T)
        Omega_new = P.dot(graph.Omega).dot(P.T)
    
        graph_new = graph_from_parameters(Lambda_new, Omega_new, prime=graph.prime)
        graph_list.append(graph_new)

    return graph_list


def vizualize_graph(graph):
    g = Digraph('G', filename="graph.gv")
    added_pairs = []
    
    for variable in graph.variables:
        g.node(variable.name)
        
        for parent in variable.parents:
            g.edge(parent.name, variable.name)
        for sibling in variable.siblings:
            if not (variable, sibling) in added_pairs:
                g.edge(variable.name, sibling.name, dir='both', style="dotted")
                added_pairs.append((sibling, variable))
                added_pairs.append((variable, sibling))
    return g


def graph_from_parameters(Lambda, Omega, prime):
    alphabet = ["", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"]
    num_variables = Lambda.shape[0]
    nodes = []
    omitted_nodes = []

    for i in range(num_variables):
        node_name = alphabet[i // 26] + alphabet[(i % 26) + 1]
        nodes.append(node_name)

    edges = []
    omitted_count = 0
    
    for row_id in range(num_variables):
        for col_id in range(num_variables):
            if Lambda[row_id, col_id]:
                edges.append((nodes[row_id], nodes[col_id]))

            if col_id > row_id:
                if Omega[row_id, col_id]:
                    omitted_count += 1
                    omitted_nodes.append(f"o{omitted_count}")
                    edges.append((f"o{omitted_count}", nodes[row_id]))
                    edges.append((f"o{omitted_count}", nodes[col_id]))

    graph = DMG(nodes, edges, omitted_nodes, prime=prime)

    return graph