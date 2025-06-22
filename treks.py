def get_half_treks(nodeA, nodeB):
    
    half_treks = []
    paths = [[(nodeA, "s")]]

    if nodeA == nodeB:
        half_treks.append([(nodeA, "s")])
    
    for sibling in nodeA.siblings:
        if sibling == nodeB:
            half_treks.append([(nodeA, "s"), (nodeB, "b")])
        else:
            paths.append([(nodeA, "s"), (sibling, "b")])

    while paths:
        path = paths.pop(0)
        
        last_node = path[-1][0]
        
        for child in last_node.children:
            if child == nodeB:
                half_treks.append(path + [(child, "r")])

            elif not (child, "r") in path:
                paths.append(path + [(child, "r")])

    return half_treks


def get_treks(nodeA, nodeB):
    treks = []
    paths = [[(nodeA, "s")]]

    # Right arrow paths:
    while paths:
        path = paths.pop(0)
        last_node = path[-1][0]

        for half_trek in get_half_treks(last_node, nodeB):
            if len(half_trek) > 1:
                treks.append(path + half_trek[1:]) 

        for parent in last_node.parents:
            if parent == nodeB:
                treks.append(path + [(parent, "l")])

            elif not (parent, "l") in path:
                paths.append(path + [(parent, "l")])

    return treks


def tr(graph, node):
    tr = []

    for var in graph.variables:
        treks = get_treks(node, var)

        if treks:
            tr.append(var)

    if not node in tr:
        tr.append(node)
        
    return tr
        
        
def htr(graph, node):
    htr = []

    for var in graph.variables:
        treks = get_half_treks(node, var)

        if treks:
            htr.append(var)
    
    if not node in htr:
        htr.append(node)

    return htr 