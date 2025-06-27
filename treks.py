from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from Variable import Variable
    from DMG import DMG 

TrekList = list[list[tuple["Variable", str]]]

def get_half_treks(nodeA: "Variable", nodeB: "Variable"):
    """
        Compute a list of tuples representing half-treks from node A to node B.
        Each tuple consists of a node and a label indicating the type of edge used to reach that node.
            - "s" = Start node
            - "r" = An edge directed to the right
            - "b" = Bidirected edge

        Args:
            nodeA (Variable): Start node.
            nodeB (Variable): End node.

        Returns:
            half_treks (TrekList): A list of half-treks from node A to node B.
    """
    half_treks = []                                                     # Initialize result list
    paths = [[(nodeA, "s")]]                                            # Add start node to the queue

    # If start and end node are the same, add trivial half-trek from the node to itself
    if nodeA == nodeB:                                  
        half_treks.append([(nodeA, "s")])

    # Start optionally with a bi-directed edge
    for sibling in nodeA.siblings:
        # If siblings is end node, append path to half_treks
        if sibling == nodeB:
            half_treks.append([(nodeA, "s"), (nodeB, "b")])             
        else:
            paths.append([(nodeA, "s"), (sibling, "b")])

    # Explore all paths
    while paths:
        path = paths.pop(0)                                             # Dequeue next path
        last_node = path[-1][0]                                         # Identify in which node we currently are
        
        for child in last_node.children:    
            new_step = (child, "r")
            # If child is end node, append path to half_treks                            
            if child == nodeB:
                half_treks.append(path + [new_step])                

            # If an edge hasn't been visited yet, append path to queue
            elif new_step not in path:
                paths.append(path + [new_step])

    return half_treks


def get_treks(nodeA, nodeB):
    """
        Compute a list of tuples representing treks from node A to node B.
        Each tuple consists of a node and a label indicating the type of edge used to reach that node.
            - "s" = Start node
            - "r" = An edge directed to the right
            - "l" = An edge directed to the left
            - "b" = Bidirected edge

        Args:
            nodeA (Variable): Start node.
            nodeB (Variable): End node.

        Returns:
            treks (TrekList): A list of treks from node A to node B.
    """
    treks = []                                                          # Initialize result list
    paths = [[(nodeA, "s")]]                                            # Add start node to the queue

    # If start and end node are the same, add trivial trek from the node to itself
    if nodeA == nodeB:                                  
        treks.append([(nodeA, "s")])

    while paths:
        path = paths.pop(0)                                             # Dequeue next path
        last_node = path[-1][0]                                         # Identify in which node we currently are

        # Add all half-treks to the list of treks
        for half_trek in get_half_treks(last_node, nodeB):
            if len(half_trek) > 1:
                treks.append(path + half_trek[1:]) 

        # If parent is end node, append path to list of treks      
        for parent in last_node.parents:
            new_step = (parent, "l")
            if parent == nodeB:
                treks.append(path + [new_step])

            # If an edge hasn't been visited yet, append step to queue
            elif new_step not in path:
                paths.append(path + [new_step])

    return treks


def tr(graph: "DMG", node: "Variable"):
    """
       Compute a list of nodes that are trek reachable from the input node

        Args:
            graph (DMG): The DMG to search in.
            node (Variable): The input node

        Returns:
            tr (list[Variable]): A list of nodes trek reachable from input node
    """
    tr = []                                                 # Initialize list of trek-reachable nodes

    for var in graph.variables:
        treks = get_treks(node, var)                        # Obtain all treks from input node to var

        if treks:
            tr.append(var)                                  # If at least one trek exists, var is trek-reachable

    # Ensure the input node is included (A node is always trek-reachable from itself)
    if node not in tr:
        tr.append(node)
        
    return tr
        
        
def htr(graph: "DMG", node: "Variable"):
    """
       Compute a list of nodes that are half-trek reachable from the input node

        Args:
            graph (DMG): The DMG to search in.
            node (Variable): The input node

        Returns:
            htr (list[Variable]): A list of nodes half-trek reachable from input node
    """
    htr = []                                                  # Initialize list of trek-reachable nodes

    for var in graph.variables:
        treks = get_half_treks(node, var)                    # Obtain all treks from input node to var

        if treks:                                           # If at least one trek exists, var is trek-reachable
            htr.append(var)
    
    # Ensure the input node is included (A node is always trek-reachable from itself)
    if node not in htr:
        htr.append(node)

    return htr 