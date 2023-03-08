from numpy import delete, zeros
from arbre import BinaryTree

class Upgma() :
    
    def __init__(self) :
        """_summary_
        """
        self.mini = 100000
        self.line_min = 0
        self.col_min = 0
        self.matrice_temp = []
    
    def find_closer_upgma(self, matrice_distance):
        """find the two clades to merge first by finding
        the minimal distance in the given matrix.

        Args:
            matrix_distance (list): matrix of the distances

        Returns:
            int: the minimum of the all the score
        """
        # Go through the mat to find the minimum :
        for x in range(0, len(matrice_distance)):
            for y in range(0, len(matrice_distance)):
                # Not comparing a seq with itself !
                if x != y:
                    # Temp_min = a value at axes x and y
                    temp_min = matrice_distance[x, y]
                    # New temp min if less than previous one
                    if temp_min < self.mini and temp_min != 0 :
                        # udpate the mini and its localisation in the matrix
                        self.mini, self.line_min, self.col_min = temp_min, x+1, y+1

    def update_upgma(self, matrice_dist) :
        """_summary_

        Args:
            matrice_dist (_type_): _description_
        """
        # Initialisation of the matrix
        size_matrix = len(matrice_dist)
        self.matrice_temp = zeros((size_matrix+1, size_matrix+1))
        # Copy the dist_matrix
        self.matrice_temp[1:,1:] = matrice_dist
        # For each group -> mean of the two values we merge.
        for i in range(1, size_matrix+1) :
            self.matrice_temp[0, i] = (self.matrice_temp[self.line_min+1, i] + self.matrice_temp[self.col_min+1, i])/2
        # Delete the line of one of the two groups we merged.
        self.matrice_temp = delete(self.matrice_temp,[self.line_min+1, self.col_min+1],axis=0)
    

    def tree_with_upgma(self, dist_matrix, clades_names) :
        """creates a tree in newick format thank to a matrix distance using upgma
        This function call the find_cloder_upgma function

        Args:
            self.matrice_distance (matrix): matrix of the distances
        
        Returns:
            tree in newick format
        """
        clades = [BinaryTree([name,[],[]]) for name in clades_names]
        length = len(clades_names)
        for i in range(length-1):
            self.find_closer_upgma(dist_matrix)
            print(dist_matrix)
            clades[self.line_min].branch_length = self.mini/2 - clades[self.line_min].depth
            clades[self.col_min].branch_length = self.mini/2 - clades[self.col_min].depth
            new_clade = clades[self.line_min].join(clades[self.col_min])
            new_clade.depth = self.mini/2
            dist_matrix = self.update_upgma(dist_matrix)
        return clades[0]
       
       
# def wpgma(distance,names):
# 	clades = [bt.BinaryTree([name,[],[]]) for name in names]
# 	n = len(names)
# 	for c in range(n-1):
# 		mindist,i,j = find_wpgma(distance)
# 		print(distance)
# 		clades[i].branch_length = mindist/2 - clades[i].depth
# 		clades[j].branch_length = mindist/2 - clades[j].depth
# 		u = clades[i].join(clades[j])
# 		# Il manque mettre a jour les valeurs dans u 
# 		u.depth = mindist/2
# 		#u.left.branch_length = u.depth - u.left.depth
# 		#u.right.branch_length = u.depth - u.right.depth
# 		distance = update_wpgma(distance,i,j)
# 		# mettre a jour clades
# 		clades = [u] + clades[:i]+clades[i+1:j]+clades[j+1:]
# 	return clades[0]

if __name__ == "__main__":
    upgma = Upgma()


# def find_wpgma(distance):
# 	n = len(distance)
# 	mindist = distance[0,1]
# 	mini,minj = 0,1
# 	for i in range(n-1):
# 		for j in range(i+1,n):
# 			if distance[i,j] < mindist:
# 				mindist,mini,minj = distance[i,j],i,j
# 	return mindist,mini,minj



# def wpgma(distance,names):
# 	clades = [bt.BinaryTree([name,[],[]]) for name in names]
# 	n = len(names)
# 	for c in range(n-1):
# 		mindist,i,j = find_wpgma(distance)
# 		print(distance)
# 		clades[i].branch_length = mindist/2 - clades[i].depth
# 		clades[j].branch_length = mindist/2 - clades[j].depth
# 		u = clades[i].join(clades[j])
# 		# Il manque mettre a jour les valeurs dans u 
# 		u.depth = mindist/2
# 		#u.left.branch_length = u.depth - u.left.depth
# 		#u.right.branch_length = u.depth - u.right.depth
# 		distance = update_wpgma(distance,i,j)
# 		# mettre a jour clades
# 		clades = [u] + clades[:i]+clades[i+1:j]+clades[j+1:]
# 	return clades[0]

# if __name__ == "__main__":
# 	# On definit la matrice de distance
# 	d = np.array([[0,5,9,9,8],
# 					[5,0,10,10,9],
# 					[9,10,0,8,7],
# 					[9,10,8,0,3],
# 					[8,9,7,3,0]])
# 	# On va donner un nom a chaque espece
# 	names = ["a","b","c","d","e"]

# 	t = wpgma(d,names)
# 	print(t.newick())