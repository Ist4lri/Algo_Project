class BinayrTree():
    """To make the Binary Tree"""

    def __init__(self) -> None:
        self.root = None

    def tree_basic(self) -> str:
        """Print the entire Tree (and not at a node)"""
        return self.root.display_node()

    def tree_verticallu(self) -> str:
        """Print the entire Tree vertically"""
        return self.root.display_vertically_node()

    def tree_depth(self) -> int:
        """Return the total depth of the tree"""
        return self.root.get_max_depth()

    def tree_node(self) -> int:
        """Return the total number of nodes in tree"""
        return self.root.node_number()

    def tree_difficult(self) -> str:
        """Trying the hardway print tree"""
        return self.root.display_hard_mode(self.root.display_hard_node())


class Node():
    """To make the node of the Tree"""

    def __init__(self, value) -> None:
        self.right = None
        self.left = None
        self.depth = 0
        self.value = value

    def link_node(self, left=None, right=None) -> None:
        """Add child node to parent node"""
        self.right = right
        self.left = left
        if self.left:
            left.depth = self.depth + 1
            left.update_children_depth()
        if self.right:
            right.depth = self.depth + 1
            right.update_children_depth()

    def update_children_depth(self) -> None:
        """Update the depth of the child node"""
        if self.left:
            self.left.depth = self.depth + 1
            self.left.update_children_depth()
        if self.right:
            self.right.depth = self.depth + 1
            self.right.update_children_depth()

    def __str__(self) -> str:
        """To have a V/D form"""
        return str(self.value) + "/" + str(self.depth)

    def node_is_leaf(self) -> bool:
        """No child node ?"""
        return not (self.right or self. left)

    def get_max_depth(self, depth_max=0) -> int:
        """Return the max depth of node"""
        if self.node_is_leaf():
            if self.depth > depth_max:
                return self.depth
        if self.right:
            depth_max = self.right.get_max_depth(depth_max)
        if self.left:
            depth_max = self.left.get_max_depth(depth_max)
        return depth_max

    def display_node(self) -> str:
        """Simplye display of nodes"""
        tabulation = ""
        for _ in range(0, self.depth):
            tabulation += "\t"
        result = tabulation + str(self) + "\n"
        if self.right:
            result += self.right.display_node()
        if self.left:
            result += self.left.display_node()
        return result

    def display_vertically_node(self) -> str:
        """Display nodes vertically"""
        max_depth = self.get_max_depth()
        tabulation = ""
        for _ in range(0, max_depth):
            tabulation += "\t"
        result = tabulation + str(self)
        if self.right and self.left:
            result += "\n" + str(self.left.display_vertically_node()) +\
                str(self.right.display_vertically_node())
        if self.left and not self.right:
            result += "\n" + str(self.left.display_vertically_node())
        if self.right and not self.left:
            result += "\n" + tabulation + \
                str(self.right.display_vertically_node())
        return result

    def node_number(self) -> int:
        """Count the number of node"""
        count = 1
        if self.node_is_leaf():
            return 1
        if self.right and not self.left:
            count += self.right.node_number()
        if self.left and not self.right:
            count += self.left.node_number()
        if self.left and self.right:
            count += self.left.node_number()
            count += self.right.node_number()
        return count

    def display_hard_mode(self, dict) -> str:
        max_depth = self.get_max_depth()
        tabulation = ""
        for _ in range(0, max_depth):
            tabulation += "\t"
        result = tabulation + \
            str(list(dict.keys())[0]) + " - " + \
            str(list(dict.values())[0]) + "\n"
        count = 0
        for key, value in dict.items():
            count += 1
            if value < list(dict.values())[0]:
                result += key + " - " + str(value)
            elif value > list(dict.values())[0]:
                result += tabulation + key + " - " + str(value)
            if (count % 3) == 0:
                result += "\n"

        return result
