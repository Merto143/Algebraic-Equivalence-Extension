import galois
import numpy as np

from Identifiability import HTC_identifiable, EID_identifiable, TSID_identifiable, EID_TSID_identifiable
from treks import htr
from Variable import Variable

class DMG():
    def __init__(self, nodes, edges, omitted_nodes = [], prime = 2**31-1):
        """
        Initialize a Directed Mixed Graph (DMG) object.

        Args:
            nodes (list): The list of observed variables (nodes).
            edges (list): A list of edges (directed or bidirected) defining the graph structure.
            omitted_nodes (list, optional): Latent variable list.
            prime (int, optional): Prime number to define the finite field GF(p) for algebraic computations.
        """
        self.nodes = nodes
        self.omitted_nodes = omitted_nodes
        self.all_nodes = nodes + omitted_nodes
        self.edges = edges
        self.prime = prime
        self.GF = galois.GF(self.prime)               # Define the finite field GF(p) used for algebraic operations

        # Track number of nodes in different categories
        self.N_nodes = len(nodes)
        self.N_omitted = len(self.omitted_nodes)
        self.N_total = len(self.all_nodes)
        
        # Store identifiability results 
        self.HTC_identifiable = None
        self.EID_identifiable = None
        self.TSID_identifiable = None
        self.EID_TSID_identifiable = None
        
        # Mappings from Variable objects and names to indices
        self.var_to_int = {}
        self.str_to_var = {}

         # Initialize zero matrices for model parameters over GF(p)
        self.Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))
        self.extended_Lambda = self.GF.Zeros((self.N_total, self.N_total))
        self.Omega = self.GF.Zeros((self.N_nodes, self.N_nodes))
        self.Sigma = self.GF.Zeros((self.N_nodes, self.N_nodes))

        # Build graph structure
        self.create_variables()     # Create Variable objects for all nodes
        self.set_parent_child()     # Establish parent/child relationships
        self.set_siblings()         # Establish sibling relationships 
        self.set_descendants()      # Obtain each variables descendant sets
        self.reset_id_grid()        # Initialize matrix for tracking which edges have been identified 

    
    def create_variables(self):
        """
        Create Variable objects for all observed and omitted nodes in the graph.

        This method:
            - Instantiates Variable objects for each observed node with a unique index.
            - Instantiates Variable objects for omitted nodes.
            - Populates internal mappings and stores the full list of Variable objects.

        Returns:
            list: The list of Variable objects corresponding to observed nodes.
        """
        # Create a variable for each of the nodes
        self.variables = []
        self.omitted_variables = []

        # Create and index Variable objects for observed nodes
        for index, node in enumerate(self.nodes):
            variable = Variable(node, index=index)
            self.variables.append(variable)
            self.var_to_int[variable] = index
            self.str_to_var[node] = variable

         # Create Variable objects for omitted nodes
        for o_node in self.omitted_nodes:
            variable = Variable(o_node, True)
            self.omitted_variables.append(variable)

        # Combine observed and omitted variables into a unified list
        self.all_variables = self.variables + self.omitted_variables

        return self.variables

    
    def set_parent_child(self):
        """
        Establish parent-child relationships between variables based on the graph's directed edges.

        For each edge (A → B) in the edge list:
            - A becomes a parent of B.
            - B becomes a child of A.

        Returns:
            list: The list of observed Variable objects (self.variables).
        """
        # Iterate over all pairs of variables (including omitted ones)
        for row_var in self.all_variables:
            for col_var in self.all_variables:
                for edge in self.edges:
                    # If edge is directed: row_var → col_var establish parent-child relationship
                    if row_var.name == edge[0] and col_var.name == edge[1]:
                        col_var.add_parent(row_var)
                        row_var.add_child(col_var)

        return self.variables


    def set_siblings(self):
        """
        Establish sibling relationships between observed variables.

        Returns:
            list: The list of observed Variable objects (self.variables).
        """

        # Iterate over each unique pair of observed variables
        for var1_id, var1 in enumerate(self.variables):
            for var2 in self.variables[var1_id + 1:]:
                # Check for shared omitted parents
                for o_parent in var1.o_parents:
                    if o_parent in var2.o_parents:
                        var1.add_sibling(var2)      # Mark var1 and var2 as siblings 
                        var2.add_sibling(var1)

        return self.variables

    
    def set_descendants(self):
        """
        Compute and assign the set of descendants for each observed variable.

        Returns:
            list: The list of observed Variable objects (self.variables).
        """
        for var in self.variables:
            var.find_descendants()      # Populate var.descendants

        return self.variables

    
    def reset_id_grid(self):
        """
        Reset the edge ID grid for directed edges in the graph.

        The grid is a matrix where entry (i, j) represents the identifiability status
        of the edge from variable i to variable j:
            - 0  → no edge
            - -1 → edge exists and is currently unidentifiable
            - 1  → edge has been identified

        This method:
            - Initializes a square matrix with zeros.
            - Sets -1 for all directed edges that exist in the graph.

        Returns:
            np.ndarray: The identifiability grid (id_grid).
        """
        self.id_grid = np.zeros((self.N_nodes, self.N_nodes))   # Start with all zeros
        
         # For each directed edge var → child, mark as unidentifiable (-1)
        for var_id, var in enumerate(self.variables):
            for child in var.children:
                child_id = self.var_to_int[child]
                self.id_grid[var_id, child_id] = -1

        return self.id_grid
                

    def generate_parameters(self):  
        """
        Generate the structural parameters of the model:
            - Lambda: directed edge coefficients
            - Omega: bidirected edge variances/covariances
            - Sigma: implied covariance matrix
        """
        self.gen_Lambda()
        self.gen_Omega()
        self.gen_Sigma()

    
    def gen_Lambda(self):
        """
        Generate the Lambda matrix, which contains coefficients for directed edges.

        For each directed edge (parent → child), a random non-zero value in GF(p) is assigned.
        This is done in the extended Lambda matrix (includes omitted nodes), and then
        the top-left submatrix is extracted for use with observed variables only.

        Returns:
            np.ndarray: The Lambda matrix (observed portion).
        """
        for var_id, var in enumerate(self.all_variables):
            for child in var.children:
                child_id = self.var_to_int[child]
                lambd = self.GF(np.random.randint(1, self.prime))       # Random value in GF(p)
                self.extended_Lambda[var_id][child_id] = lambd          # Assign to extended matrix

         # Extract observed part of Lambda for use in Σ calculation
        self.Lambda = self.extended_Lambda[:self.N_nodes, :self.N_nodes]

        return self.Lambda

    
    def gen_Omega(self):
        """
        Generate the Omega matrix, which contains variances and covariances for bidirected edges.

        - Diagonal entries represent variances and are assigned random non-zero values in GF(p).
        - Off-diagonal entries are symmetric and represent covariances between siblings.

        Returns:
            np.ndarray: The Omega matrix.
        """
        for var_id, var in enumerate(self.variables):
            # Assign random non-zero variance to the diagonal
            variance = self.GF(np.random.randint(1, self.prime))
            self.Omega[var_id][var_id] = variance
            
            # Assign symmetric covariances for each sibling pair
            for sibling in var.siblings:
                cov = self.GF(np.random.randint(1, self.prime))
                sib_id = self.var_to_int[sibling]
                self.Omega[var_id][sib_id] = cov
                self.Omega[sib_id][var_id] = cov

        return self.Omega

    
    def gen_Sigma(self):
        """
        Compute the implied covariance matrix from Lambda and Omega.

        Uses the formula:
            Sigma = (I - Lambda)^(-T) * Omega * (I - Lambda)^(-1)

        If Lambda or Omega has not been generated yet, they are generated first.

        Returns:
            np.ndarray: The implied covariance matrix Σ.
        """

        if not self.Lambda.any():
            self.gen_Lambda()
        if not self.Omega.any():
            self.gen_Omega()

        # Compute (I - Lambda) and its inverse over GF(p)
        ILambd = self.GF.Identity(self.N_nodes) - self.Lambda
        ILambd_inv = np.linalg.inv(ILambd)

        # Compute Sigma = (I - Lambda)^(-T) * Omega * (I - Lambda)^(-1)
        self.Sigma = np.matmul(ILambd_inv.T, np.matmul(self.Omega, ILambd_inv))
        
        return self.Sigma
        

    def Lambda_HTC(self, Sigma: np.ndarray, order: list, Lambda: bool = None):
        """
        Compute the Lambda matrix using the Half-Trek Criterion (HTC) identifiability method.

        This implementation follows the alorithm from the Appendix of Foygel et al. (2012):

        Given an HTC-identifiable graph and an ordering over the nodes (with their identifying sets Yv),
        this method constructs the Lambda matrix column by column by solving systems of linear equations.

        Args:
            Sigma (np.ndarray): The covariance matrix of the model.
            order (list): The HTC identification order as a list of (v, Yv) pairs.
                        - v: a node in the graph.
                        - Yv: an HTC-identifying set for v.
            Lambda (np.ndarray, optional): If provided, start from a pre-initialized Lambda matrix.

        Returns:
            Tuple[np.ndarray, dict]:
                - Lambda: The completed Lambda matrix 
                - det_multiplier: Dictionary mapping each node v to the determinant of its system matrix A.
                                This is used to compute Lambda_tilde later.
        """
        # Initialize Lambda matrix if not given
        if Lambda == None:
            Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes)) 

        det_multiplier = {}                                             # Store determinant multipliers for use in Lambda_tilde
        I = self.GF.Identity(self.N_nodes)                        
        
        # Process each (v, Yv) pair in the HTC order
        for item in order:
            var, Yv = item
            var_id = self.var_to_int[var]

            # Initialize system A x = b to solve for Lambda values into var
            A = self.GF.Zeros((len(var.parents), len(var.parents)))
            b = self.GF.Zeros(len(var.parents))
            
            # Compute (I - Lambda)ᵀ × Sigma 
            ILambdT_Sig = np.matmul((I - Lambda).T, Sigma)
            htr_var = htr(self, var)                                    # Get half-trek-reachable nodes from var
            
            # Fill in matrix A and vector b row by row
            for Yvar_index, Yvar in enumerate(Yv):
                Yvar_id = self.var_to_int[Yvar]
                
                for parent_index, parent in enumerate(var.parents):
                    parent_id = self.var_to_int[parent]
                    
                    if Yvar in htr_var:      
                        # If Yvar is trek-reachable from var
                        A[Yvar_index, parent_index] = ILambdT_Sig[Yvar_id, parent_id]
                    else:
                        A[Yvar_index, parent_index] = Sigma[Yvar_id, parent_id]

                # Fill b similarly, depending on whether Yvar in htr(var)
                if Yvar in htr_var:
                    b[Yvar_index] = ILambdT_Sig[Yvar_id, var_id]
                    
                else:
                    b[Yvar_index] = Sigma[Yvar_id, var_id]
                    
            # Solve A × Lambda_col = b for Lambda entries corresponding to var
            if var.parents:
                det = np.linalg.det(A)
                
                # If we find a zero-denominator, return an empty Lambda an det_multiplier
                if det == 0:
                    Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))
                    det_multiplier = {}

                    return Lambda, det_multiplier

                # Solve the linear system and assign Lambda values
                Lambda_new = np.matmul(np.linalg.inv(A), b)
    
                for parent_index, parent in enumerate(var.parents):
                    parent_id = self.var_to_int[parent]
                    Lambda[parent_id, var_id] = Lambda_new[parent_index]

                # Store determinant used for normalizing Lambda_tilde later
                det_multiplier[var] = det

        return Lambda, det_multiplier
        

    def Lambda_EID(self, Sigma: np.ndarray, order: list, Lambda: bool = None):
        """
        Compute the Lambda matrix using the Edgewise Identification (EID) method.

        This function implements the algorithm proposed by Weihs et al. (2018).

        The method solves, for each identified edge subset W → v, a system of equations derived
        from the covariance matrix Sigma and the structure of the graph, incorporating known values
        from previously identified edges.

        Args:
            Sigma (np.ndarray): The covariance matrix of the model.
            order (list): The EID identification order — a list of triples (v, W, Y), where:
                - v: A node whose incoming edges from W are to be identified.
                - W: A subset of v's parents to be identified.
                - Y: A set of observed variables used as an EID-identifying set for W → v.
            Lambda (np.ndarray, optional): An optional partially filled Lambda matrix.

        Returns:
            Tuple[np.ndarray, dict]:
                - Lambda: The estimated Lambda matrix after identification.
                - det_multiplier: Dictionary mapping v to the determinant of the linear system
                                used to solve for its Lambda column.
        """
        if not np.any(Lambda):
            Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))         # Initialize if not provided
            
        det_multiplier = {}                                              # Store determinant multipliers for normalization in Lambda_tilde

        # Process each (v, W, Y) triple in the EID order
        for v, W, Y in order:    
            A = self.GF.Zeros((len(W), len(W)))
            b = self.GF.Zeros(len(W))
            v_id = self.var_to_int[v]

            # Construct the linear system A * Λ_Wv = b
            for y_index, y in enumerate(Y):
                y_id = self.var_to_int[y]

                H = []                                                   # H = known parents of y (edges into y already identified)
                for parent in y.parents:
                    parent_id = self.var_to_int[parent]

                    if not Lambda[parent_id, y_id] == 0:
                        H.append(parent)

                # Fill row y_index of matrix A
                for w_index, w in enumerate(W):
                    w_id = self.var_to_int[w]
                    h_sum_1 = self.GF(0)
                    
                    for h in H:
                        h_id = self.var_to_int[h]
                        h_sum_1 += Sigma[h_id, w_id]*Lambda[h_id, y_id]

                    A[y_index, w_index] = Sigma[y_id, w_id] - h_sum_1

                # Compute RHS entry for b[y_index]
                p_sum = self.GF(0)
                for p in (set(v.parents) - set(W)):
                    p_id = self.var_to_int[p]

                    h_sum_2 = self.GF(0)
                    for h in H:
                        h_id = self.var_to_int[h]
                        h_sum_2 += Sigma[h_id, p_id]*Lambda[h_id, y_id]

                    if (Sigma[y_id, p_id] - h_sum_2):   
                        p_sum += (Sigma[y_id, p_id] - h_sum_2)*Lambda[p_id, v_id]

                h_sum_3 = self.GF(0)
                for h in H:
                    h_id = self.var_to_int[h] 
                    h_sum_3 += Sigma[v_id, h_id]*Lambda[h_id, y_id]
      
                b[y_index] = Sigma[y_id, v_id] - p_sum - h_sum_3
            
            # Obtain denominator of the rational function for Lambda
            det = np.linalg.det(A)

            # If we find a zero-denominator, return an empty Lambda an det_multiplier
            if det == 0:
                Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))
                det_multiplier = {}

                return Lambda, det_multiplier
            
            else:
                # Solve for the new Lambda values
                A_inv = np.linalg.inv(A)
                            
                new_lambda = np.matmul(A_inv, b)

                for w_index, w in enumerate(W):
                    w_id = self.var_to_int[w]
                    Lambda[w_id, v_id] = new_lambda[w_index]

            # Track determinant multiplier for finding Lambda_tilde later
            if det_multiplier.get(v, None):
                det_multiplier[v] = det_multiplier[v] * det

            else:
                det_multiplier[v] = det

        return Lambda, det_multiplier
        

    def Lambda_TSID(self, Sigma: np.ndarray, order: list, Lambda: bool = None):
        """
        Compute the Lambda matrix using the Trek-Separation Identification (TSID) method.

         This function implements the algorithm proposed by Weihs et al. (2018).

        For each edge w₀ → v that is identified via trek-separation, this method uses determinant
        ratios involving subsets of source (S) and sink (T) variables to solve for the corresponding
        Lambda coefficient.

        Args:
            Sigma (np.ndarray): The observed covariance matrix of the model.
            order (list): The TSID identification order — a list of (w₀, v, S, T) tuples, where:
                - w₀: A parent of v to be identified,
                - v: The child node,
                - S: A set of source nodes,
                - T: A set of sink nodes.
            Lambda (np.ndarray, optional): A pre-initialized Lambda matrix. If not provided, zeros are used.

        Returns:
            Tuple[np.ndarray, dict]:
                - Lambda: Updated Lambda matrix with newly identified entries from TSID.
                - det_multiplier: Dictionary mapping each node v to the product of determinants used for normalization.
        """
        
        if not np.any(Lambda):
            Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))    # Initialize if needed

        det_multiplier = {}                                         # Save determinants for each node

        # Process each TSID-identifiable edge (w0 → v)
        for item in order:
            w0, var, S_comb, T_comb = item
            w0_id = self.var_to_int[w0]
            var_id = self.var_to_int[var]
            
            # Convert variable objects to matrix indices
            S_ids = [self.var_to_int[key] for key in S_comb]
            T_ids = [self.var_to_int[key] for key in T_comb]

            # Compute det1: det(Σ[S, T ∪ {v}])
            col_ids = T_ids + [var_id]
            det1 = np.linalg.det(Sigma[S_ids, :][:, col_ids])

            # Compute det2: det(Σ[S, T ∪ {w₀}])
            col_ids = T_ids + [w0_id]
            det2 = np.linalg.det(Sigma[S_ids, :][:, col_ids])

            lambda_sum = self.GF(0)         # Sum over already identified parents ≠ w0
            
            # Subtract contributions from previously identified parents of v
            for parent in var.parents:
                parent_id = self.var_to_int[parent]
                if self.id_grid[parent_id, var_id] == 1 and parent != w0:
                    lamb = Lambda[parent_id, var_id]
                    col_ids = T_ids + [parent_id]
                    det3 = np.linalg.det(Sigma[S_ids, :][:, col_ids])

                    lambda_sum = lambda_sum + lamb*det3

            # If we find a zeo denominator
            if det2 == 0:
                Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))
                det_multiplier = {}

                return Lambda, det_multiplier

            # Solve for the Lambda entry 
            Lambda[w0_id, var_id] = (det1 - lambda_sum)/det2

            # Update determinant multiplier for normalization in Lambda_tilde
            if det_multiplier.get(var, None):
                det_multiplier[var] = det_multiplier[var] * det2
    
            else:
                det_multiplier[var] = det2
            
        return Lambda, det_multiplier
            
    
    def Lambda_EID_TSID(self, Sigma: np.ndarray, order: list):
        """
    Compute the Lambda matrix using a combined EID + TSID identification strategy.

    This method executes both Edgewise Identification (EID) and Trek-Separation Identification (TSID)
    steps in a user-defined order (as provided in `order`), accumulating Lambda entries and tracking
    the determinant multipliers for normalization.

    Each item in the input `order` list is a tuple:
        - ("EID", EID_order): where EID_order is a list of (v, W, Y) triples.
        - ("TSID", TSID_order): where TSID_order is a list of (w0, v, S, T) tuples.

    Args:
        Sigma (np.ndarray): The observed covariance matrix of the graph.
        order (list): A sequence of identification steps, alternating between EID and TSID,
                      each represented as a ("METHOD", order_data) pair.

    Returns:
        Tuple[np.ndarray, dict]:
            - Lambda: The estimated Lambda matrix, filled by both EID and TSID methods.
            - det_multiplier: Dictionary mapping each node v to the product of determinant
                              scalings encountered during identification.
    """
        Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))        # Initialize Lambda with zeros
        det_multiplier = {}                                         # Track determinants used in Lambda_tilde computation

        # Process each ("EID" or "TSID", suborder) entry in the combined order
        for item in order:
            if item[0] == "EID":
                EID_order = item[1]
                Lambda, EID_mult = self.Lambda_EID(Sigma, EID_order, Lambda=Lambda)

                # Update determinant multipliers
                for var in EID_mult:
                    if det_multiplier.get(var, None):
                        det_multiplier[var] = det_multiplier[var] * EID_mult[var]

                    else:
                        det_multiplier[var] = EID_mult[var]

                # If EID step failed (Lambda is all zeros), return early
                if not np.any(Lambda):
                    Lambda_tilde = Lambda
            
                    return Lambda_tilde, det_multiplier
  
            elif item[0] == "TSID":
                TSID_order = item[1]
                Lambda, TSID_mult = self.Lambda_TSID(Sigma, TSID_order, Lambda=Lambda)

                # Update determinant multipliers
                for var in TSID_mult:
                    if det_multiplier.get(var, None):
                        det_multiplier[var] = det_multiplier[var] * TSID_mult[var]
    
                    else:
                        det_multiplier[var] = TSID_mult[var]
                
                # If TSID step failed (Lambda is all zeros), return early
                if not np.any(Lambda):
                    Lambda_tilde = Lambda
            
                    return Lambda_tilde, det_multiplier

        return Lambda, det_multiplier

    
    def Lambda_tilde(self, Sigma: np.ndarray, method: str):
        """
        Compute the transformed Lambda matrix using an identification method.

        The resulting Lambda_tilde is used in containment and identifiability checks.

        Args:
            Sigma (np.ndarray): The covariance matrix of the graph.
            method (str): The identification method to use. Options are:
                        "HTC", "EID", "TSID", "EID_TSID".

        Returns:
            np.ndarray: The matrix Lambda_tilde used in computing Omega_tilde.
        """
        Lambda_tilde = self.GF.Identity(self.N_nodes)       # Start with identity matrix 
        
        # === HTC method ===
        if method == "HTC":
            # Run HTC identifiability check if not already cached
            if self.HTC_identifiable == None:
                self.HTC_identifiable, self.HTC_order = HTC_identifiable(self)

            # Raise an exception if HTC fails
            if self.HTC_identifiable == False:
                raise Exception("Graph not HTC identifiable")

            # Compute Lambda and its determinant multipliers
            Lambda, det_multiplier = self.Lambda_HTC(Sigma, self.HTC_order)

        # === EID method ===
        elif method == "EID":
            if self.EID_identifiable == None:
                self.EID_identifiable, self.EID_order = EID_identifiable(self)
    
            if self.EID_identifiable == False:
                raise Exception("Graph not EID identifiable")
            
            Lambda, det_multiplier = self.Lambda_EID(Sigma, self.EID_order)

        # === TSID method ===
        elif method == "TSID":
            if self.TSID_identifiable == None:
                self.TSID_identifiable, self.TSID_order = TSID_identifiable(self)
    
            if self.TSID_identifiable == False:
                raise Exception("Graph not TSID identifiable")

            Lambda, det_multiplier = self.Lambda_TSID(Sigma, self.TSID_order)

        # === Combined EID + TSID method ===
        elif method == "EID_TSID":
            if self.EID_TSID_identifiable == None:
                self.EID_TSID_identifiable, self.EID_TSID_order = EID_TSID_identifiable(self)
            
            if self.EID_TSID_identifiable == False:
                raise Exception("Graph not EID+TSID identifiable")
                
            Lambda, det_multiplier = self.Lambda_EID_TSID(Sigma, self.EID_TSID_order)

        else:
            raise Exception("Choose a valid identificaton method")

        # For each variable v, adjust the corresponding column in Lambda tilde
        for v in det_multiplier:
            v_id = self.var_to_int[v]
            Lambda_tilde[:, v_id] = (Lambda_tilde[:, v_id] - Lambda[:, v_id]) * det_multiplier[v]

        # If a zero denominator has been found, return zero matrix
        if not np.any(Lambda):
            Lambda_tilde = Lambda
                        
        return Lambda_tilde

        
    def Omega_tilde(self, Sigma, Lambda_tilde):
        """
        Compute the implied Omega matrix Omega_tilde from the given Sigma and Lambda_tilde.

        Args:
            Sigma (np.ndarray): The covariance matrix Σ of the model being tested.
            Lambda_tilde (np.ndarray): The estimated Lambda matrix Λ̃ under the current method.

        Returns:
            np.ndarray: The implied bidirected edge matrix Ω̃.
        """
        Omega_tilde = np.matmul(Lambda_tilde.T, np.matmul(Sigma, Lambda_tilde))
        
        return Omega_tilde