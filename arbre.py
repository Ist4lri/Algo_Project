class BinaryTree:
    def __init__(self, l):
        self.newick_dict = {}

        if len(l) == 0:
            self.val = None
            self.nLeaf = 0
            self.depth = 0
            self.branch_length = 0
            self.left = None
            self.right = None
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
                self.depth = 1 + max(self.left.depth, self.right.depth)

    def is_terminal(self):
        """Check if we are on a leaf or not

        Returns:
            Boolean: Are we on,a leaf ?
        """
        return (self.val is None)

    def join(self, rightTree):
        """For Joining the Tree

        Args:
            rightTree (List): List of BinaryTree

        Returns:
            Class: BinaryTree
        """
        tmp = BinaryTree([])
        tmp.val = ""
        tmp.branch_length = 0
        tmp.left = self
        tmp.right = rightTree
        tmp.nLeaf = self.nLeaf + rightTree.nLeaf
        tmp.depth = 0
        return tmp

    def newick(self):
        """To make a Newick Format

        Returns:
            String: The BinaryTree in Newick Format
        """
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

    def parse_newick(self, newick):
        """Fonction pour parser le newick en un arbre sous forme de dictionnaire

        Args:
            newick (string): The newick format of tree

        Returns:
            dictionnary: The Newick in dictionnary
        """
        newick = newick.replace("(", "{").replace(")", "}")
        newick = newick.split(";")[0]
        newick_list = newick.split(",")
        self.newick_dict = {}
        for node in newick_list:
            node_name = node.split(":")[0].replace(
                "(", "").replace(")", "").replace("{", "").replace("}", "")
            node_dist = float(node.split(":")[1].replace(",", ".").rstrip("}"))
            if "{" in node_name:
                clade1, clade2 = node_name.split("{")[1].split("}")
                clade1_name, clade1_dist = clade1.split(":")
                clade2_name, clade2_dist = clade2.split(":")
                self.newick_dict[node_name] = {"clades": [clade1_name, clade2_name], "distances": [
                    float(clade1_dist), float(clade2_dist)], "dist": node_dist}
            else:
                self.newick_dict[node_name] = node_dist
        return self.newick_dict
