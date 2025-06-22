import networkx as nx

def max_flow(graph, v, A):
    g = nx.DiGraph()

    g.add_edge("s", "t", capacity=0)

    for var in A:            
        g.add_edge("s", f"l_{var.name}_in")
        g.add_edge(f"l_{var.name}_out", f"r_{var.name}_in")
    
        for sibling in var.siblings:
            g.add_edge(f"l_{var.name}_out", f"r_{sibling.name}_in")

    for var in graph.variables:
        g.add_edge(f"l_{var.name}_in", f"l_{var.name}_out", capacity=1)
        g.add_edge(f"r_{var.name}_in", f"r_{var.name}_out", capacity=1)
        
        for child in var.children: 
            g.add_edge(f"r_{var.name}_out", f"r_{child.name}_in")

    for parent in v.parents:
        g.add_edge(f"r_{parent.name}_out", "t")

    flow_value, flow_dict = nx.maximum_flow(g, 's', 't')
    
    return flow_value, flow_dict


def max_flow_TSID(graph, S, T, remove_edges=[]):
    g = nx.DiGraph()

    g.add_edge("s", "t", capacity=0)
    
    for var in S:            
        g.add_edge("s", f"l_{var.name}_in")

    for var in T:
        g.add_edge(f"r_{var.name}_out", "t")
    
    for var in graph.variables:
        g.add_edge(f"l_{var.name}_in", f"l_{var.name}_out", capacity=1)
        g.add_edge(f"r_{var.name}_in", f"r_{var.name}_out", capacity=1)

        g.add_edge(f"l_{var.name}_out", f"r_{var.name}_in")

        for child in var.children:
            g.add_edge(f"l_{child.name}_out", f"l_{var.name}_in")
            g.add_edge(f"r_{var.name}_out", f"r_{child.name}_in")
        
        for sibling in var.siblings:
            g.add_edge(f"l_{var.name}_out", f"r_{sibling.name}_in")
    
    for edge in remove_edges:
        parent = edge[0]
        child = edge[1]

        g.remove_edge(f"r_{parent.name}_out", f"r_{child.name}_in")

    flow_value, flow_dict = nx.maximum_flow(g, 's', 't')

    
    return flow_value, flow_dict
    

