from numpy import delete, zeros
from arbre import BinaryTree


class Upgma():

    def __init__(self):
        """_summary_
        """
        self.mini = 0
        self.line_min = 0
        self.col_min = 0

    def find_closer_upgma(self, matrice_distance):
        """find the two clades to merge first by finding
        the minimal distance in the given matrix.

        Args:
            matrix_distance (list): matrix of the distances

        Returns:
            int: the minimum of the all the score
        """
        self.mini = 10000
        self.line_min = 0
        self.col_min = 0
        # Go through the mat to find the minimum :
        for x in range(0, len(matrice_distance)):
            for y in range(0, len(matrice_distance)):
                # Temp_min = a value at axes x and y
                temp_min = matrice_distance[x, y]
                # New temp min if less than previous one
                if temp_min < self.mini and temp_min != 0:
                    # udpate the mini and its localisation in the matrix
                    self.mini, self.line_min, self.col_min = temp_min, x, y

    def update_upgma(self, matrice_dist):
        """_summary_

        Args:
            matrice_dist (_type_): _description_
        """
        # Initialisation of the matrix
        size_matrix = len(matrice_dist)
        matrice_temp = zeros((size_matrix+1, size_matrix+1))
        # Copy the dist_matrix
        matrice_temp[1:, 1:] = matrice_dist
        # For each group -> mean of the two values we merge.
        for i in range(1, size_matrix+1):
            matrice_temp[i, 0] = (matrice_temp[self.line_min+1, i] +
                                  matrice_temp[self.col_min+1, i])/2
            matrice_temp[0, i] = matrice_temp[i, 0]
        # Delete the line of one of the two groups we merged.
        matrice_temp = delete(
            matrice_temp, [self.line_min+1, self.col_min+1], axis=0)
        matrice_temp = delete(
            matrice_temp, [self.line_min+1, self.col_min+1], axis=1)
        return matrice_temp

    def tree_with_upgma(self, dist_matrix, clades_names):
        """creates a tree in newick format thank to a matrix distance using upgma
        This function call the find_cloder_upgma function

        Args:
            self.matrice_distance (matrix): matrix of the distances
            clades_names (list) : names of the species

        Returns:
            tree in newick format, list of clades in order
        """
        # On crée un arbre binaire
        clades = [BinaryTree([name, [], []]) for name in clades_names]
        length = len(clades_names)
        # On trouve les distances les plus courtes entre les séquences des espèces de la liste
        for i in range(length-1):
            self.find_closer_upgma(dist_matrix)
            # On calcule la distance (et localisation dans la matrice) entre les deux espèces les plus proches
            clades[self.line_min].branch_length = self.mini / \
                2 - clades[self.line_min].depth
            clades[self.col_min].branch_length = self.mini / \
                2 - clades[self.col_min].depth
            new_clade = clades[self.line_min].join(clades[self.col_min])
            new_clade.depth = self.mini/2
            # Tant qu'il reste des taxons
            if len(dist_matrix) > 1:
                dist_matrix = self.update_upgma(dist_matrix)
                clades = [new_clade] + clades[:self.line_min] + \
                    clades[self.line_min+1:self.col_min] + \
                    clades[self.col_min+1:]
        return clades[0].parse_newick(clades[0].newick()), clades[0].newick()
