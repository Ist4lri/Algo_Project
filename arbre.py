class BinaryTree:
    def __init__(self,l):
        if len(l) == 0:
            self.val = None
            self.nLeaf = 0
            self.depth = 0
            self.branch_length = 0
            self.left=None
            self.right=None
        else:
            self.val = l[0]
            self.branch_length = 0
            self.left = BinaryTree(l[1])
            self.right = BinaryTree(l[2])
            if self.left.is_terminal() and self.right.is_terminal():
                self.nLeaf = 1
                self.depth = 0
            else:
                self.nLeaf = self.left.nLeaf + self.right.nLeaf 
                self.depth = 1 + max(self.left.depth,self.right.depth)

    def is_terminal(self):
        return (self.val is None)

    def join(self,rightTree):
        tmp = BinaryTree([])
        tmp.val = ""
        tmp.branch_length = 0
        tmp.left = self
        tmp.right = rightTree
        tmp.nLeaf = self.nLeaf + rightTree.nLeaf
        tmp.depth = 0
        return tmp

    def newick(self):
        if self.left.is_terminal() and self.right.is_terminal():
            return f"{self.val}:{self.branch_length:.2f}"
        return f'({self.left.newick()},{self.right.newick()}):{self.branch_length:.2f}'
    
    def __str__(self):
        if self.is_terminal():
            return ''
        else: 
            s = str(self.val)
            if self.left is not None:
                s = s + str(self.left)
            if self.right is not None:
                s = s + str(self.right)
            return s

    
if __name__ == "__main__":
    # On teste la crÃ©ation d'un arbre
    bt = BinaryTree([1, [2, [4,[],[5,[],[]] ] ,[] ], [3,[],[]]])
    print(bt)