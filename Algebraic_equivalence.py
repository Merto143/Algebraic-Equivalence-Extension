import numpy as np

def model_contained(G, G_prime, method):
    """
    Check whether model G is contained in model G' using the specified identification method.

    This function uses the parameter matrices (Lambda_tilde and Omega_tilde) derived from G and G' to determine
    whether the covariance structure of G can be represented by G'.

    Args:
        G (DMG): The input model whose structure is being tested for containment.
        G_prime (DMG): The comparison model.
        method (str): The identification method to use (e.g., "EID", "EID_TSID", "HTC").

    Returns:
        bool or str:
            - True if G is contained in G'.
            - False if G is not contained in G'.
            - "INC" if the method finds a zero-denominator.
    """
    contained = True                                                # Assume containment unless proven otherwise
    G.generate_parameters()                                         # Generate covariance matrix from Gs parameters
    Lambda_tilde = G_prime.Lambda_tilde(G.Sigma, method)            # Compute Lambda tilde from G' based on Gs covariance matrix

    # If all Lambda entries are zero, a zero-denominator has been found
    if np.all(Lambda_tilde == 0):
        return "INC"

    # Compute Omega tilde from G' using Lambda tilde and Gs covariance matrix
    Omega_tilde = G_prime.Omega_tilde(G.Sigma, Lambda_tilde)

    # Check that Omega Tilde only has nonzero entries for bidirected edges in G'
    for var in G_prime.variables:
        for var_2 in G_prime.variables:
            if not var_2 in var.siblings and var != var_2:
                var_id = G_prime.var_to_int[var]
                var_2id = G_prime.var_to_int[var_2]

                if not Omega_tilde[var_id, var_2id] == 0:
                    contained = False                                # Nonzero entry where no bidirected edge exists → not contained
                     
    return contained

def algebraic_equivalence(graph_a, graph_b, method):
    """
    Check whether two graphs are algebraically equivalent using a specified identification method.

    Two graphs are algebraically equivalent if each model is algebraically contained in the other.

    Args:
        graph_a (DMG): The first graph.
        graph_b (DMG): The second graph.
        method (str): The identification method to use (e.g., "EID", "EID_TSID", "HTC").

    Returns:
        bool or str:
            - True if the graphs are algebraically equivalent.
            - False if they are not equivalent.
            - "INC" if the result is inconclusive due to a zero denominator.
    """
    a_in_b = model_contained(graph_a, graph_b, method)      # Check if graph_a is contained in graph_b
    b_in_a = model_contained(graph_b, graph_a, method)      # Check if graph_b is contained in graph_a

     # If either comparison is inconclusive, and the other is inconclusive or True, return "INC"
    if (a_in_b == "INC" and b_in_a == "INC") or (a_in_b == "INC" and b_in_a == True) or (a_in_b == True and b_in_a == "INC"):
        return "INC"

    # If both containment checks returned True, graphs are algebraically equivalent
    elif a_in_b == True and b_in_a == True:
        equivalence = True

    # Otherwise, they are not equivalent
    else:
        equivalence = False
                                
    return equivalence

