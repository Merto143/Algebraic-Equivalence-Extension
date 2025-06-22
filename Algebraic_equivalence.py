import numpy as np

def model_contained(G, G_prime, method):
    contained = True
    G.generate_parameters()
    Lambda_tilde = G_prime.Lambda_tilde(G.Sigma, method)

    if np.all(Lambda_tilde == 0):
        return "INC"

    Omega_tilde = G_prime.Omega_tilde(G.Sigma, Lambda_tilde)

    for var in G_prime.variables:
        for var_2 in G_prime.variables:
            if not var_2 in var.siblings and var != var_2:
                var_id = G_prime.var_to_int[var]
                var_2id = G_prime.var_to_int[var_2]

                if not Omega_tilde[var_id, var_2id] == 0:
                    contained = False
                     
    return contained

def algebraic_equivalence(graph_a, graph_b, method):
    a_in_b = model_contained(graph_a, graph_b, method)
    b_in_a = model_contained(graph_b, graph_a, method)

    if (a_in_b == "INC" and b_in_a == "INC") or (a_in_b == "INC" and b_in_a == True) or (a_in_b == True and b_in_a == "INC"):
        return "INC"
                
    elif a_in_b == True and b_in_a == True:
        equivalence = True

    else:
        equivalence = False
                                
    return equivalence

