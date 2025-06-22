class Variable():
    def __init__(self, name, omitted=False, index=-1):
        self.name = name
        self.all_parents = []
        self.o_parents = []
        self.parents = []
        self.children = []
        self.siblings = []
        self.omitted = omitted
        self.Yv = []
        self.index = index
        
    def add_parent(self, parent):
        self.all_parents.append(parent)

        if parent.omitted == False:
            self.parents.append(parent)
        else:
            self.o_parents.append(parent)

        return 

    
    def add_child(self, child):
        self.children.append(child)

        return self.children

        
    def add_sibling(self, node):
        self.siblings.append(node)

        return self.siblings
        
    
    def add_Yv(self, var):
        self.Yv.append(var)

        return self.Yv

    def add_descendents(self):
        self.descendents = []
        children = self.children
        
        change = True
        while change:
            change = False
            
            for child in children:
                if not child in self.descendents:
                    self.descendents.append(child)
                    children = children + child.children
                    change = True

    def __lt__(self, other):
        return self.index < other.index
                
    def __repr__(self):
        return self.name