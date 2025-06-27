class Variable():
    """
    A class representing a variable

    Attributes:
        name (str): The name of the variable.
        omitted (bool): Whether this is an omitted variable.
        index (int): The natural order of variable (-1 if the order doesn't matter).

        parents (List[Variable]): Observed parents.
        o_parents (List[Variable]): Unobserved parents.
        all_parents (List[Variable]): All parents.
        children (List[Variable]): Children of this variable.
        siblings (List[Variable]): Siblings of this variable.
        descendants (List[Variable]): All descendants of this variable.

    Methods:
        add_parent(parent): Adds another variable to the list of parents.
        add_child(child): Adds another variable to the list of children.
        add_sibling(sibling): Adds another variable to the list of siblings.
        find_descendants(): Obtain a list containing all descendants.
    """
    def __init__(self, name: str, omitted: bool = False, index: int = -1):
        """
        Initialize a new variable

        Args:
            name (str): The name of the variable.
            omitted (bool): Whether this is an omitted variable.
            index (int): The natural order of variable (-1 if the order doesn't matter).

        """
        self.name = name
        self.omitted = omitted
        self.index = index

        self.parents = [] 
        self.o_parents = [] 
        self.all_parents = [] 
        
        self.children = [] 
        self.siblings = [] 

    def add_parent(self, parent: "Variable"):
        """
        Add a variable to the appropriate parent list.

        Args:
            parent (Variable): The variable to add to list of parents.
        """
        self.all_parents.append(parent)    # Add to the list of all parents

        if parent.omitted:
            self.o_parents.append(parent)  # Add to the list of unobserved parents
        else:
            self.parents.append(parent)    # Add to the list of observed parents 

    
    def add_child(self, child: "Variable"):
        """
        Add a variable to the list of children.

        Args:
            child (Variable): The variable to add to list of children.
        """
        self.children.append(child)

        
    def add_sibling(self, sibling: "Variable"):
        """
        Add a variable to the list of children.

        Args:
            child (Variable): The variable to add to list of children.
        """
        self.siblings.append(sibling)
        

    def find_descendants(self):
        """
        Obtain a list of descendants.

        """ 
        self.descendants = []
        visited = set()                             # Track variables already visited
        queue = list(self.children)                 # Queue current children
        
        while queue:  
            child = queue.pop(0)
            if child not in visited:
                self.descendants.append(child)      # Add child to the list of descendants
                visited.add(child)                  # Mark variable as visited variables
                queue.extend(child.children)        # Add grandchildren to the queue


    def __lt__(self, other):
        return self.index < other.index
                
    def __repr__(self):
        return self.name