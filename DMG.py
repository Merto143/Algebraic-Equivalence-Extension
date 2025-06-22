import galois
import numpy as np

from Identifiability import HTC_identifiable, EID_identifiable, TSID_identifiable, EID_TSID_identifiable
from treks import htr
from Variable import Variable

class DMG():
    def __init__(self, nodes, edges, omitted_nodes = [], prime = 2**31-1):
        self.nodes = nodes
        self.omitted_nodes = omitted_nodes
        self.all_nodes = nodes + omitted_nodes
        self.edges = edges
        self.prime = prime
        self.GF = galois.GF(self.prime)

        self.N_nodes = len(nodes)
        self.N_omitted = len(self.omitted_nodes)
        self.N_total = len(self.all_nodes)
        
        self.HTC_identifiable = None
        self.EID_identifiable = None
        self.TSID_identifiable = None
        self.EID_TSID_identifiable = None
        
        self.var_to_int = {}
        self.str_to_var = {}

        # Initialize empty parameters
        self.Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))
        self.extended_Lambda = self.GF.Zeros((self.N_total, self.N_total))
        self.Omega = self.GF.Zeros((self.N_nodes, self.N_nodes))
        self.Sigma = self.GF.Zeros((self.N_nodes, self.N_nodes))

        # Create graph structure
        self.create_variables()
        self.set_parent_child()
        self.set_siblings()
        self.set_descendents()
        self.reset_id_grid()

    
    def create_variables(self):
        # Create a variable for each of the nodes
        self.variables = []
        self.omitted_variables = []

        
        for index, node in enumerate(self.nodes):
            # instantiate observed variables
            variable = Variable(node, index=index)
            self.variables.append(variable)
            self.var_to_int[variable] = index
            self.str_to_var[node] = variable

        for o_node in self.omitted_nodes:
            # instantiate omitted variables
            variable = Variable(o_node, True)
            self.omitted_variables.append(variable)

        self.all_variables = self.variables + self.omitted_variables

        return self.variables

    
    def set_parent_child(self):
        for row_id, row_var in enumerate(self.all_variables):
            for col_id, col_var in enumerate(self.all_variables):
                for edge in self.edges:
                    # if row variable is a parent of the column variable
                    if row_var.name == edge[0] and col_var.name == edge[1]:
                        # set parent-child relation and set the lambda value
                        col_var.add_parent(row_var)
                        row_var.add_child(col_var)

        return self.variables


    def set_siblings(self):
        # go trough each pair of variables
        for var1_id, var1 in enumerate(self.variables):
            for var2 in self.variables[var1_id + 1:]:
                for o_parent in var1.o_parents:
                    # for each variable, check if they have the same omitted parents
                    if o_parent in var2.o_parents:
                        # set sibling status for both variables
                        var1.add_sibling(var2)
                        var2.add_sibling(var1)

        return self.variables

    
    def set_descendents(self):
        for var in self.variables:
            var.add_descendents()

        return self.variables

    
    def reset_id_grid(self):
        self.id_grid = np.zeros((self.N_nodes, self.N_nodes))
        
        for var_id, var in enumerate(self.variables):
            for child in var.children:
                child_id = self.var_to_int[child]
                self.id_grid[var_id, child_id] = -1

        return self.id_grid
                

    def generate_parameters(self):  
        self.gen_Lambda()
        self.gen_Omega()
        self.gen_Sigma()

    
    def gen_Lambda(self):
        for var_id, var in enumerate(self.all_variables):
            for child in var.children:
                child_id = self.var_to_int[child]
                lambd = self.GF(np.random.randint(1, self.prime))
                self.extended_Lambda[var_id][child_id] = lambd

        self.Lambda = self.extended_Lambda[:self.N_nodes, :self.N_nodes]

        return self.Lambda

    
    def gen_Omega(self):
        for var_id, var in enumerate(self.variables):
            variance = self.GF(np.random.randint(1, self.prime))
            self.Omega[var_id][var_id] = variance

            for sibling in var.siblings:
                cov = self.GF(np.random.randint(1, self.prime))
                sib_id = self.var_to_int[sibling]
                self.Omega[var_id][sib_id] = cov
                self.Omega[sib_id][var_id] = cov

        return self.Omega

    
    def gen_Sigma(self):
        if not self.Lambda.any():
            self.gen_Lambda()
        if not self.Omega.any():
            self.gen_Omega()

        ILambd = self.GF.Identity(self.N_nodes) - self.Lambda
        ILambd_inv = np.linalg.inv(ILambd)
        self.Sigma = np.matmul(ILambd_inv.T, np.matmul(self.Omega, ILambd_inv))
        
        return self.Sigma
        

    def Lambda_HTC(self, Sigma, order, Lambda = None):
        if Lambda == None:
            Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes)) 

        det_multiplier = {}
        I = self.GF.Identity(self.N_nodes)
        
        for item in order:
            var, Yv = item
            var_id = self.var_to_int[var]
            A = self.GF.Zeros((len(var.parents), len(var.parents)))
            b = self.GF.Zeros(len(var.parents))
            
            ILambdT_Sig = np.matmul((I - Lambda).T, Sigma)
            htr_var = htr(self, var)
            
            for Yvar_index, Yvar in enumerate(Yv):
                Yvar_id = self.var_to_int[Yvar]
                
                for parent_index, parent in enumerate(var.parents):
                    parent_id = self.var_to_int[parent]
                
                    if Yvar in htr_var:      
                        A[Yvar_index, parent_index] = ILambdT_Sig[Yvar_id, parent_id]
                
                    else:
                        A[Yvar_index, parent_index] = Sigma[Yvar_id, parent_id]

                if Yvar in htr_var:
                    b[Yvar_index] = ILambdT_Sig[Yvar_id, var_id]
                    
                else:
                    b[Yvar_index] = Sigma[Yvar_id, var_id]
                    
            
            if var.parents:
                det = np.linalg.det(A)

                if det == 0:
                    Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))
                    det_multiplier = {}

                    return Lambda, det_multiplier

                Lambda_new = np.matmul(np.linalg.inv(A), b)
    
                for parent_index, parent in enumerate(var.parents):
                    parent_id = self.var_to_int[parent]
                    Lambda[parent_id, var_id] = Lambda_new[parent_index]
    
                det_multiplier[var] = det

        return Lambda, det_multiplier
        

    def Lambda_EID(self, Sigma, order, Lambda = None):    
        if not np.any(Lambda):
            Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))
            
        det_multiplier = {}

        for v, W, Y in order:    
            # Calculate the A matrix
            A = self.GF.Zeros((len(W), len(W)))
            b = self.GF.Zeros(len(W))
            v_id = self.var_to_int[v]

            for y_index, y in enumerate(Y):
                y_id = self.var_to_int[y]

                H = []
                for parent in y.parents:
                    parent_id = self.var_to_int[parent]

                    if not Lambda[parent_id, y_id] == 0:
                        H.append(parent)

                for w_index, w in enumerate(W):
                    w_id = self.var_to_int[w]
                    h_sum_1 = self.GF(0)

                    for h in H:
                        h_id = self.var_to_int[h]
                        h_sum_1 += Sigma[h_id, w_id]*Lambda[h_id, y_id]

                    A[y_index, w_index] = Sigma[y_id, w_id] - h_sum_1

                # calculate rhs of the equation:
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

            det = np.linalg.det(A)

            if det == 0:
                for h in H:
                    h_id = self.var_to_int[h]
                Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))
                det_multiplier = {}
                
                return Lambda, det_multiplier
                
            A_inv = np.linalg.inv(A)
                        
            new_lambda = np.matmul(A_inv, b)

            for w_index, w in enumerate(W):
                w_id = self.var_to_int[w]
                lambd = new_lambda[w_index]
                Lambda[w_id, v_id] = new_lambda[w_index]

            if det_multiplier.get(v, None):
                det_multiplier[v] = det_multiplier[v] * det

            else:
                det_multiplier[v] = det

        return Lambda, det_multiplier
        

    def Lambda_TSID(self, Sigma, order, Lambda = None):
        if not np.any(Lambda):
            Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))

        det_multiplier = {}

        for item in order:
            w0, var, S_comb, T_comb = item
            w0_id = self.var_to_int[w0]
            var_id = self.var_to_int[var]
            
            S_ids = [self.var_to_int[key] for key in S_comb]
            T_ids = [self.var_to_int[key] for key in T_comb]
        
            col_ids = T_ids + [var_id]
            det1 = np.linalg.det(Sigma[S_ids, :][:, col_ids])
        
            col_ids = T_ids + [w0_id]
            det2 = np.linalg.det(Sigma[S_ids, :][:, col_ids])

            lambda_sum = self.GF(0)
            
            for parent in var.parents:
                parent_id = self.var_to_int[parent]
                if self.id_grid[parent_id, var_id] == 1 and parent != w0:
                    lamb = Lambda[parent_id, var_id]
                    col_ids = T_ids + [parent_id]
                    det3 = np.linalg.det(Sigma[S_ids, :][:, col_ids])

                    lambda_sum = lambda_sum + lamb*det3

            if det2 == 0:
                Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes))
                det_multiplier = {}

                return Lambda, det_multiplier

            Lambda[w0_id, var_id] = (det1 - lambda_sum)/det2

            if det_multiplier.get(var, None):
                det_multiplier[var] = det_multiplier[var] * det2
    
            else:
                det_multiplier[var] = det2
            
        return Lambda, det_multiplier
            
    
    def Lambda_EID_TSID(self, Sigma, order):
        Lambda = self.GF.Zeros((self.N_nodes, self.N_nodes)) 
        det_multiplier = {}

        for item in order:
            if item[0] == "EID":
                EID_order = item[1]
                Lambda, EID_mult = self.Lambda_EID(Sigma, EID_order, Lambda=Lambda)

                for var in EID_mult:
                    if det_multiplier.get(var, None):
                        det_multiplier[var] = det_multiplier[var] * EID_mult[var]
    
                    else:
                        det_multiplier[var] = EID_mult[var]

                if not np.any(Lambda):
                    Lambda_tilde = Lambda
            
                    return Lambda_tilde, det_multiplier
  
            elif item[0] == "TSID":
                TSID_order = item[1]
                Lambda, TSID_mult = self.Lambda_TSID(Sigma, TSID_order, Lambda=Lambda)

                for var in TSID_mult:
                    if det_multiplier.get(var, None):
                        det_multiplier[var] = det_multiplier[var] * TSID_mult[var]
    
                    else:
                        det_multiplier[var] = TSID_mult[var]
                
                if not np.any(Lambda):
                    Lambda_tilde = Lambda
            
                    return Lambda_tilde, det_multiplier

        return Lambda, det_multiplier

    
    def Lambda_tilde(self, Sigma, method):
        Lambda_tilde = self.GF.Identity(self.N_nodes)
        
        # Check if the graph is HTC identifiable
        if method == "HTC":
            if self.HTC_identifiable == None:
                self.HTC_identifiable, self.HTC_order = HTC_identifiable(self)

            # Raise exception for non-HTC identifiable graphs
            if self.HTC_identifiable == False:
                raise Exception("Graph not HTC identifiable")

            Lambda, det_multiplier = self.Lambda_HTC(Sigma, self.HTC_order)


        elif method == "EID":
            # Check if the graph is EID identifiable
            if self.EID_identifiable == None:
                self.EID_identifiable, self.EID_order = EID_identifiable(self)
    
            # Raise exception for non-EID identifiable graphs
            if self.EID_identifiable == False:
                raise Exception("Graph not EID identifiable")
            
            Lambda, det_multiplier = self.Lambda_EID(Sigma, self.EID_order)


        elif method == "TSID":
            # Check if the graph is TSID identifiable
            if self.TSID_identifiable == None:
                self.TSID_identifiable, self.TSID_order = TSID_identifiable(self)
    
            # Raise exception for non-TSID identifiable graphs
            if self.TSID_identifiable == False:
                raise Exception("Graph not TSID identifiable")

            Lambda, det_multiplier = self.Lambda_TSID(Sigma, self.TSID_order)


        elif method == "EID_TSID":
            # Check if the graph is EID_TSID identifiable
            if self.EID_TSID_identifiable == None:
                self.EID_TSID_identifiable, self.EID_TSID_order = EID_TSID_identifiable(self)
            
            # Raise exception for non-EID_TSID identifiable graphs
            if self.EID_TSID_identifiable == False:
                raise Exception("Graph not EID+TSID identifiable")
                
            Lambda, det_multiplier = self.Lambda_EID_TSID(Sigma, self.EID_TSID_order)

        else:
            raise Exception("Choose a valid identificaton method")

        
        for v in det_multiplier:
            v_id = self.var_to_int[v]
            Lambda_tilde[:, v_id] = (Lambda_tilde[:, v_id] - Lambda[:, v_id]) * det_multiplier[v]

        
        if not np.any(Lambda):
            Lambda_tilde = Lambda
                        
        return Lambda_tilde

        
    def Omega_tilde(self, Sigma, Lambda_tilde):
        Omega_tilde = np.matmul(Lambda_tilde.T, np.matmul(Sigma, Lambda_tilde))
        
        return Omega_tilde