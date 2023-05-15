from constants import *
from numpy import array, delete, zeros, ones, transpose, unravel_index, argmin
from numpy import sum as np_sum
from numpy import max as np_max
from upgma import Upgma
from arbre import BinaryTree


class Alignement():

    def __init__(self):
        self.upgma = Upgma()
        self.clades_names = []  # Liste des noms de séquences (en seqX)
        self.nj_clades_group = []
        self.matrice_distance = []  # Matrice de distance
        self.dir_matrix = []  # Matrice des directions (Après 2e NW)
        self.dict_tree = {}  # Dico après premier UPGMA arbre {'seqX' : distance}
        self.conserved_dict = {}
        # Dico avec les séquences pour alignement {'seqX' :  AGTCGAT...}
        self.dict_sequence_tree = {}
        self.alignement1 = []  # Premier alignement
        self.alignement2 = []  # Deuxième alignement si nécessaire
        self.conserved_alignement = []
        self.dist_mat_conserved= [] # matrice des distances entre les sequences conservées
        self.freq_mat1 = []  # matrice de fréquences pour l'alignement1
        self.freq_mat2 = []  # matrice de fréquences pour l'alignement2


    def needleman_wunsch(self, seq1, seq2, matrix, gap=-1) -> int:
        """
        Algorithme de Needleman-Wunsch pour aligner deux séquences

        Arguments:
        seq1 -- première séquence (str)
        seq2 -- deuxième séquence (str)
        gap -- score de trou (int, default=-1)

        Retourne:
        score -- score d'alignement (int)
        """
        # Initialisation de la matrice de scores
        n = len(seq1)
        m = len(seq2)
        score_matrix = [[0] * (m+1) for i in range(n+1)]
        # Initialisation de la première ligne et colonne de la matrice
        for i in range(1, n+1):
            score_matrix[i][0] = gap * i
        for j in range(1, m+1):
            score_matrix[0][j] = gap * j
        # Remplissage de la matrice de scores
        for i in range(1, n+1):
            for j in range(1, m+1):
                indice_lig = DAA[seq1[i-1]]
                indice_col = DAA[seq2[j-1]]
                match_mismatch_score = matrix[indice_lig, indice_col]
                score_matrix[i][j] = max(
                    score_matrix[i-1][j-1] +
                    match_mismatch_score,  # Correspondance
                    score_matrix[i-1][j] + gap,  # Trou dans seq1
                    score_matrix[i][j-1] + gap  # Trou dans seq2
                )
        return score_matrix[n][m]

    def get_all_max_score(self, fasta_dict):
        """get the max score of all sequence

        Args:
            filename (str): path of the file

        Returns:
            list: array of array with all score
        """
        counter_col = 0
        counter_lig = 0
        index_dict = len(fasta_dict)
        # Put only zeros in a mat of the good size
        self.matrice_distance = array([[0] * (index_dict)
                                       for _ in range(index_dict+1)])
        # shorter names :
        self.clades_names = fasta_dict.keys()
        # For each pair of seq :
        for header in fasta_dict.keys():
            counter_col = 0
            counter_lig += 1
            seq_one = fasta_dict[header]
            for key, seq_two in fasta_dict.items():
                if header != key:
                    # find the score between both seq and put it in the distance mat :
                    self.matrice_distance[counter_lig, counter_col] = self.needleman_wunsch(
                        seq_one, seq_two, BLOSUM62)
                counter_col += 1
        # To delete the first line of zeros :
        self.matrice_distance = delete(self.matrice_distance, 0, axis=0)
        # Transform Needleman-Wunsch matrix to an Upgma matrix :
        len_mat_dist = len(self.matrice_distance)
        self.matrice_distance = self.matrice_distance - \
            np_max(self.matrice_distance)*ones((len_mat_dist, len_mat_dist))
        for line in range(len(self.matrice_distance)):
            for col in range(len(self.matrice_distance)):
                if self.matrice_distance[line, col] != 0:
                    self.matrice_distance[line, col] = - \
                        self.matrice_distance[line, col]
        #self.matrice_distance = self.upgma.transform_mat_dist(self.matrice_distance)
        # To UPGMA :
        self.upgma_to_multiple_align(self.matrice_distance, self.clades_names)

    def upgma_to_multiple_align(self, distance_matrix, clades_names):
        self.dict_tree = self.upgma.tree_with_upgma(
            distance_matrix, clades_names)
        self.multiple_alignement()

    def freq_matrix_calc(self, seq_list):
        """Cette fonction prend une liste de séquences d'acides aminés et retourne une matrice de fréquence"""
        num_seq = len(seq_list)
        seq_len = len(seq_list[0])
        freq_matrix = zeros((21, seq_len))
        # Liste des acides aminés
        amino_acids = "CSTAGPDEQNHRKMILVWYF-"
        for i in range(num_seq):
            for j in range(seq_len) :
                indice_aa_matrice = amino_acids.find(seq_list[i][j])
                freq_matrix[indice_aa_matrice][j] += 1
        freq_matrix /= num_seq
        # # Liste des acides aminés
        # amino_acids = "CSTAGPDEQNHRKMILVWYF-"
        # # Initialisation de la matrice avec des zéros
        # num_seq = len(seq_list)
        # seq_len = len(seq_list[0])
        # freq_matrix = zeros((seq_len, len(amino_acids)))
        # # Remplissage de la matrice de fréquences
        # for i in range(seq_len):
        #     for j in range(num_seq):
        #         indice_aa_dans_matrice = amino_acids.find(seq_list[j][i])
        #         freq_matrix[i][indice_aa_dans_matrice] += 1
        # freq_matrix /= num_seq
        return freq_matrix

    def needleman_wunsch_profile(self, seq1 : list, seq2 : list, matrix : list):
        """Cette fonction prend en entrée deux profils
        de séquences multiples, une matrice de score,
        et une valeur de gap (par défaut -4) et
        renvoie un nouveau profil aligné"""

        if type(seq1) is not list:
            seq1 = [seq1]
        if  type(seq2) is not list:
            seq2 = [seq2]
        
        # On trouve la séquence le profil le plus court (+ petite en ligne, + grande en colonnes)
        if len(seq1[0]) > len(seq2[0]):
            nb_col = len(seq1[0]) +1
            nb_ligne = len(seq2[0]) +1
            shorter_profile = seq2
            longer_profile = seq1
            freq_matrix1 = self.freq_matrix_calc(seq1)
            freq_matrix2 = self.freq_matrix_calc(seq2)
        else :
            nb_col = len(seq2[0]) +1
            nb_ligne = len(seq1[0]) +1
            shorter_profile = seq1
            longer_profile = seq2
            freq_matrix2 = self.freq_matrix_calc(seq1)
            freq_matrix1 = self.freq_matrix_calc(seq2)

        # définition du gap :
        gap = zeros(21)
        gap[-1] = 1
        # Initialisation de la matrice des scores et de direction
        score_matrix = zeros((nb_ligne, nb_col)) # m = nb de lignes, n = nb de colonnes
        direction_matrix = zeros((nb_ligne, nb_col))

        # On rempli la première colonne en fonction des gaps:
        for i in range(1, nb_ligne):
            # on prend le score de la case de gauche,
            # on l'additionne à la np.array([0,0,0, ...,1]) @ BLOSUM62 @ fréq
            score_matrix[i][0] = score_matrix[i-1][0] + (gap @ matrix @ freq_matrix2)[i-1]
            direction_matrix[i][0] = 1
        # On rempli la prmière ligne en fonction des gaps:
        for j in range(1, nb_col):
            # on prend le score de la case du haut,
            # on l'additionne à la fréquence d'un gap * BLOSUM62 * la fréquence observée sur la position
            score_matrix[0][j] = score_matrix[0][j-1] + (gap @ matrix @ freq_matrix1)[j-1]
            direction_matrix[0][j] = 2
        
        # On change de sens pour pouvoir faire le calcul matriciel
        freq_matrix1 = transpose(freq_matrix1)
        freq_matrix2 = transpose(freq_matrix2)
        
        # Remplissage de la matrice des scores et de direction
        for i in range(1, nb_col): # taille de la seq la plus longue (seq1)
            for j in range(1, nb_ligne): # taille de la seq la plus courte (seq2)
                # Calcul du score pour un alignement diagonal
                diag_score = score_matrix[j-1][i-1] + (freq_matrix1[i-1] @ matrix @ freq_matrix2[j-1])
                
                # Calcul du score pour un alignement horizontal
                left_score = score_matrix[j][i-1] + (freq_matrix2[j-1] @ matrix @ gap)
                #left_score = score_matrix[i-1][j] + (freq_matrix2[j-1] @ matrix @ gap)
                
                # Calcul du score pour un alignement vertical
                up_score = score_matrix[j-1][i] + (freq_matrix1[i-1] @ matrix @ gap)
                #up_score = score_matrix[i][j-1] + (freq_matrix1[i-1] @ matrix @ gap)
                
                score_matrix[j][i], direction_matrix[j][i] = max((diag_score, 0), (up_score, 1), (left_score, 2))

        # Construction du nouvel alignement
        aligned_seq1 = []
        aligned_seq2 = []

        # Pour chaque séquence de la liste/du profile
        for k in longer_profile:
            i = nb_ligne - 1
            j = nb_col - 1
            #print('\n new_aa_to_join: ',new_aa_to_join)
            new_aa_to_join = ""
            # On va chercher dans la matrice de direction
            while i > 0 or j > 0:
                direction = direction_matrix[i][j]
                # Si le score vient de la diagonale :
                if direction == 0:
                    new_letter = k[j-1]
                    new_aa_to_join += new_letter
                    i -= 1
                    j -= 1
                # Si le score vient d'en haut
                elif direction == 1:
                    new_aa_to_join += "-"
                    i -= 1
                # Si le score vient de la gauche
                elif direction == 2:
                    new_letter = k[j-1]
                    new_aa_to_join += new_letter
                    j -= 1
            aligned_seq1.append(new_aa_to_join[::-1])

        # Pour chaque séquence dans le profile
        for k in shorter_profile:
            i = nb_ligne - 1
            j = nb_col - 1
            new_aa_to_join2 = ""
            # Tant qu'il reste des AA ou nt
            while i > 0 or j > 0:
                # On regarde dans la matrice de direction
                direction = direction_matrix[i][j]
                # Si le score vient de la diagonale
                if direction == 0:
                    new_letter = k[i-1]
                    new_aa_to_join2 += new_letter
                    i -= 1
                    j -= 1
                # Si le score vient d'en haut
                elif direction == 1:
                    new_letter = k[i-1]
                    new_aa_to_join2 += new_letter
                    i -= 1
                # Si le score vient de la gauche
                elif direction == 2:
                    new_aa_to_join2 += "-"
                    j -= 1
            aligned_seq2.append(new_aa_to_join2[::-1])

        new_alignment = []

        for elem in aligned_seq1:
            new_alignment.append(elem)
        for elem in aligned_seq2:
            new_alignment.append(elem)

        return new_alignment

    def multiple_alignement(self):
        """
        Alignements multiples des séquences selon la distance trouvée avec l'arbre newick.
        Cette fonction marche !
        """
        list_seq = [] # list of clades sorted from shorter distance to longer
        list_dist = [] # list of the distances from shorter to longer
        list_boolean = []
        # On trie le dico avec en clé les séquences et en valeur la distance (les plus proches en premiers)
        for k, v in sorted(self.dict_tree.items(), key=lambda x: x[1]):
            del self.dict_tree[k]
            self.dict_tree[k] = v
        # On place les clés et valeurs dans deux listes séparées (c'est trié) :
        for key, value in self.dict_tree.items():
            list_seq.append(key)
            list_dist.append(value)
        # On change clades_names --> il est maintenant trié dans l'ordre croissant des distances
        self.clades_names = list_seq
        for _ in list_dist:
            list_boolean.append(False)
        # On sort les séquences à merge du bon dico qui les contient :

        i = 0
        while i < len(list_dist):
            if list_boolean[i] == False:
                # Si c'est la dernière séquence, on l'ajoute à l'alignement:
                if i == len(list_dist)-1:
                    espece_solo = self.dict_sequence_tree[list_seq[i]]
                    self.alignement1 = self.needleman_wunsch_profile(self.alignement1, espece_solo, BLOSUM62)
                    list_boolean[i] = True
                    return self.alignement1
                # Si les deux première distances (triées) sont identiques :
                elif list_dist[i] == list_dist[i+1]:
                    if self.alignement1 == []:
                        espece1 = self.dict_sequence_tree[list_seq[i]]
                        espece2 = self.dict_sequence_tree[list_seq[i+1]]
                        self.alignement1 = self.needleman_wunsch_profile(espece1, espece2, BLOSUM62)
                    else:
                        espece1 = self.dict_sequence_tree[list_seq[i]]
                        espece2 = self.dict_sequence_tree[list_seq[i+1]]
                        self.alignement2 = self.needleman_wunsch_profile(espece1, espece2, BLOSUM62)
                    list_boolean[i] = True
                    list_boolean[i+1] = True
                elif list_dist[i] != list_dist[i+1]:
                    self.alignement2 = [self.dict_sequence_tree[list_seq[i]]]
                    list_boolean[i] = True
                if self.alignement1 != [] and self.alignement2 != []:
                    self.alignement1 = self.needleman_wunsch_profile(self.alignement1, self.alignement2, BLOSUM62)
            i+=1
        return self.alignement1

    def conserved_position(self) :
        """
        Trouve les positions conservées dans les alignements multiples et reconstruit
        les séquences alignées en ne gardant que ces positions:
        Position conservée : 1 ou 2 aa (ou gap) différents maximum
        """
        pos_conserved = []
        temp_list = []
        # pour chaque position des séquences alignées
        for i in range(len(self.alignement1[0])):
            # pour chaque alignement on ajoute l'aa à cette position
            for align in self.alignement1:
                letter = align[i]
                temp_list.append(letter)
            # On regarde les positions conservées
            aa1 = temp_list[0]
            aa2 = None
            aa3 = None
            conserved = True
            # pour chaque aa de la liste temp
            for aa in temp_list :
                # Si aa2 (déjà deux aa différents)
                if aa2:
                    # Si juste deux positions conservées on continue si pas de nouvel aa
                    if aa in [aa1, aa2]:
                        continue
                    else: # si un troisième aa est trouvé, la position est non-conservée
                        aa3 = aa
                        conserved = False
                    # Si un seul aa pour le moment, on continue et on en ajoute un si nouveau
                if aa2 == None:
                    if aa == aa1:
                        continue
                    else:
                        aa2 = aa
            # Si un aa différent ou 2 : on ajoute la position conservée à la liste
            if conserved == True:
                pos_conserved.append(i)
            # On vide la liste temporaire
            temp_list = []

        # On reconstruit les alignements en fonctions des positions conservées
        for align in self.alignement1:
            conserved_align = ""
            for pos in pos_conserved:
                conserved_align += align[pos]
            # On extrait les mêmes alignements mais avec seulement les positions conservées (clades_names est tjr dans l'ordre)
            self.conserved_alignement.append(conserved_align)
       
                
            # comparer la position
            # s'il y en a que deux différents, on garde
            # Refaire les séquences avec juste les positions conservées

    def conserved_distance_matrix(self):
        """Created a matrix distance from the conserved sequences after multiple alignement.
        The distance between two species will be the number of differences between the conserved
        sequences of those two species."""
        # On crée un dictionnaire avec les espèces en clés et les séquences conservées en valeurs.
        for i in range(len(self.clades_names)):
            self.conserved_dict[self.clades_names[i]] = self.conserved_alignement[i]
        # Création d'un matrice de distance vide (qui ne contient que des zéros)
        conserved_dist_mat = array([[0] * len(self.conserved_dict) for _ in range(len(self.conserved_dict)+1)])
        
        # On rempli la matrice en comptant les différences entre deux séquences
        counter_col = 0
        counter_lig = 0
        for clade in self.conserved_dict.keys():
            counter_col = 0
            counter_lig += 1
            seq_one = self.conserved_dict[clade]
            for key, seq_two in self.conserved_dict.items():
                if clade != key:
                    # On compte les différences entre les deux seq et on l'ajoute à la matrice:
                    score_dist = self.count_differences_in_seq(seq_one, seq_two)
                    conserved_dist_mat[counter_lig, counter_col] = score_dist
                counter_col += 1
        # To delete the first line of zeros :
        conserved_dist_mat = delete(conserved_dist_mat, 0, axis=0)
        self.dist_mat_conserved = conserved_dist_mat
        
    def count_differences_in_seq(self, seq1, seq2):
        """Compte les différences entre deux séquences"""
        # On parcours les séquences (de même taille)
        counter = 0
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]: # si pas le même aa
                counter += 1 # la distance augmente de 1
        return counter


    def nj_matQ_creation(self):
        """Calcul de la matrice Q qui sert à ...
        """
        n = len(self.dist_mat_conserved)
        # Calcul des sommes de distances pour chaque clade
        total_dist = np_sum(self.dist_mat_conserved, axis=1)
        # Calcul de la matrice Q
        Q = zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                Q[i, j] = (n - 2) * self.dist_mat_conserved[i, j] - total_dist[i] - total_dist[j]
                Q[j, i] = Q[i, j]
        return Q, total_dist
    
    
    def nj_matrix_update(self, matrice_dist, i, j):
        """On update de la matrice de distance. Calcul des distances entre le noeud qui relie
        les taxons i et j et tous les autres clades.
        On met aussi à jour une liste de clades de l'ordre dans lequel les clades s'associent dans self.nj_clades_group
        en prenant les bons indices de self.clades_names (qui est dans l'ordre de la matrice)."""
        # Initialisation of the matrix
        size_matrix = len(matrice_dist)
        matrice_temp = zeros((size_matrix+1,size_matrix+1))
        print(matrice_temp)
        # # Copy the dist_matrix
        matrice_temp[1:,1:] = matrice_dist
        print("mat_dist :", matrice_dist)
        # Pour chaque taxon restant, on calcul sa distance avec le noeud et on l'ajoute à une somme.
        
        for k in range(1,size_matrix) :
            print("k: ", k, "i: ", i, "j: ", j)
            matrice_temp[k, 0] = (matrice_dist[k, i] + matrice_dist[k, j] - matrice_dist[i, j]) / 2
            print(matrice_dist)
            matrice_temp[0, k] = matrice_temp[k, 0]
            print(matrice_temp)
        matrice_temp = delete(matrice_temp, [i+1,j+1], axis=0)
        matrice_temp = delete(matrice_temp, [i+1,j+1], axis=1)
        print(matrice_temp)

        #Cette boucle parcourt tous les indices k de 0 à size_matrix - 1.
        # Pour chaque valeur de k différente de i et j, la distance entre le clade k
        # et le nœud qui relie les clades i et j est calculée et assignée à matrice_temp.
        # Les indices utilisés pour l'assignation sont ajustés en fonction de la présence de
        # i et j dans la matrice. L'équation utilisée pour le calcul de la distance est basée
        # sur la formule du Neighbor Joining.

        #Veuillez noter que cette partie du code semble incomplète et qu'il pourrait y avoir d'autres
        # lignes nécessaires pour effectuer la mise à jour complète de la matrice de distance et de la liste de clades.

                
        # # Delete the line of one of the two groups we merged.
        # matrice_temp = delete(matrice_temp,[self.line_min+1, self.col_min+1],axis=0)
        # matrice_temp = delete(matrice_temp,[self.line_min+1, self.col_min+1],axis=1)
        self.dist_mat_conserved = matrice_temp

        
    def nj_global(self):
        """_summary_
        """
        # On crée un arbre binaire
        tree_nj = [BinaryTree([name,[],[]]) for name in self.clades_names]
        for k in range(len(self.clades_names)-1):
            n = len(self.dist_mat_conserved)
            # Création de la matrice Q et récupérartion des distances totales pour chaque clades
            matQ, total_dist = self.nj_matQ_creation()
            # Trouver les indices (i, j) correspondant à la plus petite valeur de Q
            min_indices = unravel_index(argmin(matQ), matQ.shape)
            i, j = min_indices
            # On ne fait plus l'hypothèse de l'horloge moléculaire par rapport à upgma.
            # Donc on calcul les distance de nos deux clades les plus proches par rapport au noeud qui les relie.
            delta = (total_dist[i] - total_dist[j]) / (n - 2)
            limb_length_i = (self.dist_mat_conserved[i, j] + delta) / 2
            limb_length_j = self.dist_mat_conserved[i, j] - limb_length_i
            # On défini u le nouveau noeud
            node = tree_nj[i].join(tree_nj[j])
            # Calcul d'une nouvelle matrice de distance entre le noeud et tous les autres clades
            self.nj_matrix_update(self.dist_mat_conserved, i, j)
            print(self.dist_mat_conserved)
            # mise à jour des clades
            tree_nj = [node] + tree_nj[:i]+tree_nj[i+1:j]+tree_nj[j+1:]
        print(tree_nj[0])
            
            

        # ÇA ÇA VIENT DU SCRIPT UPGMA
        #     self.find_closer_upgma(dist_matrix)
        #     # On calcule la distance (et ocalisation dans la matrice) entre les deux espèces les plus proches
        #     clades[self.line_min].branch_length = self.mini/2 - clades[self.line_min].depth
        #     clades[self.col_min].branch_length = self.mini/2 - clades[self.col_min].depth
        #     new_clade = clades[self.line_min].join(clades[self.col_min])
        #     new_clade.depth = self.mini/2
        #     # Tant qu'il reste des taxons
        #     if len(dist_matrix) > 1:
        #         dist_matrix = self.update_upgma(dist_matrix)
        #         clades = [new_clade] + clades[:self.line_min]+clades[self.line_min+1:self.col_min]+clades[self.col_min+1:]
        # #clades = [new_clade] + clades[:self.line_min-1]+clades[self.line_min:self.col_min-1]+clades[self.col_min:]
        # return clades[0].parse_newick(clades[0].newick())
        
        # TODO
        # FAIRE L'UPDATE DU NEIGHBOR JOINING
        # METTRE EN PLACE L'UPDATE DES CLADES EN MÊME TEMPS
        # SORTIR UN ARBRE NEWICK
        # FAIRE UN MAIN
        # FAIRE UN FICHIER D'OUTPUT
        



    def neighbor_joining(self, dist_matrix, clades):
        """_summary_

        Args:
            dist_matrix (_type_): _description_
            clades (_type_): _description_

        Returns:
            _type_: _description_
        """
        n = len(dist_matrix)
        if n == 2:
            # Crée un arbre avec une seule branche reliant les deux clades
            tree = {'left': clades[0], 'right': clades[1], 'length': dist_matrix[0, 1]}
            return tree
        
        # Calcul des sommes de distances pour chaque clade
        total_dist = np_sum(dist_matrix, axis=1)
        
        # Calcul de la matrice Q
        Q = zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                Q[i, j] = (n - 2) * dist_matrix[i, j] - total_dist[i] - total_dist[j]
                Q[j, i] = Q[i, j]
        
        # Trouver les indices (i, j) correspondant à la plus petite valeur de Q
        min_indices = unravel_index(argmin(Q), Q.shape)
        i, j = min_indices
        
        # Calcul des distances entre les nouveaux clades
        delta = (total_dist[i] - total_dist[j]) / (n - 2)
        limb_length_i = (dist_matrix[i, j] + delta) / 2
        limb_length_j = (dist_matrix[i, j] - delta) / 2
        
        # Mise à jour de la matrice de distance avec les nouveaux clades
        new_dist_matrix = zeros((n-1, n-1))
        new_dist_matrix[:-1, :-1] = (dist_matrix[:-1, :-1] * (n - 2) - dist_matrix[i, :-1] - dist_matrix[j, :-1]) / (n - 2)
        new_dist_matrix[-1, :-1] = (dist_matrix[i, :-1] + dist_matrix[j, :-1] - dist_matrix[i, j]) / 2
        new_dist_matrix[:-1, -1] = new_dist_matrix[-1, :-1]
        
        # Suppression des lignes et colonnes correspondant aux anciens clades
        new_dist_matrix = delete(new_dist_matrix, (i, j), axis=0)
        new_dist_matrix = delete(new_dist_matrix, (i, j), axis=1)
        
        # Récursion pour construire l'arbre avec les nouveaux clades
        subtree_i = self.neighbor_joining(new_dist_matrix, clades[:i] + clades[i+1:j] + clades[j+1:])
        subtree_j = self.neighbor_joining(new_dist_matrix, clades[:i] + clades[i+1:j] + clades[j+1:])
        
        # Crée un nouvel arbre en reliant les sous-arbres aux nouveaux clades
        tree = {'left': subtree_i, 'right': subtree_j}
        tree['left']['length'] = limb_length_i
        tree['right']['length'] = limb_length_j
        print(tree)
        return tree






# Matrice de distance sur la base des conditions conservées.
# Dist entre deux seq = nombre de fois où ces seqs sont différentes.




if __name__ == "__main__":
    align = Alignement()
    #print(align.multiple_alignement())
    list_seq = ['MR----KM-E---------F-LFNS-VKPWDGPQYHIAPWAF-QAFMG-TVFL-G-FPLNAVL-VATYKKLRQPLNYILVNVSFGFL-CFSVF---VFVSC--GYFF-G-VCAEA-FLG-TA-GLVTGWSLAFLAFERYVICKPFGNFRFS-KHAV----VTWIGIGVSI--PPFFGWSRFPEGL-QCSCGPDWYTVGTK---YSESY-T-FLFICFIPL--LICFSY-QLLR--ALKAVAAQQQESATTQKAEREVSRMVVVMVGSFCVCYPYAAAMYVNN-R-NHG--LRLVTIPAFFSKSACIYNPIIY-FMNKQFA-HIMKMVC--G-KMDESS---SSQKTEVSTVSSVGP-',
'MR----KM-----------F-LFNS-V-PWDGPQYHIAPWAF-QAFMG-TVFL-G-FPLNAVL-VATYKKLRQPLNYILVNVSFGFL-CFSVF---VFVSC--GYFF-G-VCAEG-FLG-TA-GLVTGWSLAFLAFERYVICKPFGNFRFS-KHAV----VTWIGIGVSI--PPFFGWSRFPEGL-QCSCGPDWYTVGTK---YSESY-T-FLFICFIPL--LICFSY-QLLR--ALKAVAAQQQESATTQKAEREVSRMVVVMVGSFCVCYPYAAAMYVNN-R-NHG--LRLVTIPSFFSKSACIYNPIIY-FMNKQFAC-IMKMVC--G-KMDESTC--SSQKTEVSTVSSVGP-',
'M-------G----------F-LFNS-V-PWDGPQYHIAPWAF-QTFMGF-VF-CGT-PLNAVL-VATYKKLRQPLNYILVNVSFGFI-CFAVF---VFIS---GYFF-G-VCAEA-FLG-SA-GLVTGWSLAFLAFERFVICKPFGNFRFS-KHAV----VTWIGIGVSI--PPFFGWSRFPEGL-QCSCGPDWYTVGTK---YSE-YYT-FLFICFIPLF-LICFSY-QLL-G-ALRAVAAQQQESATTQKAEREVSRMVVMMVGSFCLCYPYAAAMYVNN-R-NHG--LRLVTIPAFFSKSACVYNPIIY-FMNKQFAC-IMEMVCR---KMDDSS---SSQKTEVSAVSSVGP-',
'M-------G----------F-LFNS-V-PWDGPQYHIAPWAF-QTFMGF-VF-CGT-PLNAVL-VATYKKLRQPLNYILVNVSLGFI-CFAVF---VFIS---GYFF-G-VCAEA-FLG-SA-GLVTGWSLAFLAFERFVICKPFGNFRFS-KHAV----VTWIGIGVSI--PPFFGWSRFPEGL-QCSCGPDWYTVGTK---YSE-YYT-FLFICFIPLF-LICFSY-QLLR--ALRAVAAQQQESATTQKAEREVSRMVVMMVGSFCLCYPYAAAMYVNN-R-NHG--LRLVTIPAFFSKSACVYNPIIY-FMNKQFAC-IMEMVCR---KMDDSS---SSQKTEVSAVSSVGP-',
'M-------G----------F-LFNS-V-PWDGPQYHIAPWAF-QTFMGF-VF--GT-PLNAVL-IATYKKLRQPLNYILVNISLGFI-CFSVF---VFIS---GYFF-G-VCAEG-FLG-SA-GLVTGWSLAFLAFERFVICKPFGNFRFS-KHSV----VTWIGIGVSI--PPFFGWSRYPEGL-QCSCGPDWYTVGTK---YSE-YYT-FLFICFIPL--LICFSY-QLL-G-ALRAVAAQQQESATTQKAEREVSRMVVMMVGSFCLCYPYAAAMYVNN-R-NHG--LRLVTIPAFFSKSSCVYNPIIY-FMNKQFAC-IMEMVCR---KMDDSS---SSQKTEVSTVSSVGP-',
'M-------G----------F-LFNS-V-PWDGPQYHIAPWAF-QAFMGF-VF--GT-PLNAVL-VATYKKLRQPLNYILVNVSLGFL-CFSVF---VFISC--GYFF-G-VCAEA-FLG-SA-GLVTGWSLAFLAFERYVICKPFGNIRFS-KHAV----VTWIGIGVSI--PPFFGWSRFPEGL-QCSCGPDWYTVGTK---YSE-HYT-FLFICFIPL--LICFSYFQLLR--TLRAVAAQQQESATTQKAEREVSHMVVVMVGSFCLCYPYAAAMYVNN-R-NHGY-LRLVTIPAFFSKSSCVYNPIIY-FMNKQFAC-ILEMVCR---KMDESS---GSQKTEVSTVSSVGP-',
'M-------G----------F-LFNS-V-PWDGPQYHLAPWAF-QAFMGF-VF--GT-PLNAVL-VATYKKLRQPLNYILVNVSLGFL-CFSVF---VFISC--GYFF-G-VCAEA-FLG-SA-GLVTGWSLAFLAFERYVICKPFGSIRFS-KHAV----VTWIGIGVSI--PPFFGWSRFPEGL-QCSCGPDWYTVGTK---YSE-YYT-FLFICFIPL--LICFSY-QLLR--TLRAVAAQQQESATTQKAEREVSHMVVVMVGSFCLCYPYAAAMYVNN-R-NHG--LRLVTIPAFFSKSSCVYNPIIY-FMNKQFAC-ILEMVCR---KMDESS---GSQKTEVSTVSSVGP-',
'M-------G----------F-LFDS-V-PWDGPQYHIAPWAF-QTFMGF-VF--GT-PLNAVL-IATYKKLRQPLNYILVNISLGFI-CISVF---VFIS---GYFF-G-VCAEA-FLG-SA-GLVTGWSLAFLAFERFVICKPFGNFRFS-KHAV----VTWIGIGVSI--PPFFGWSRYPEGL-QCSCGPDWYTVGTK---YSE-YYTGFLFICFIPL--LICFSY-QLL-G-ALRAVAAQQQESATTQKAEREVSRMVVVMVGSFCLCYPYAAAMYVNN-R-NHG--LRLVTIPAFFSKSACVYNPIVYWFMNKQFAC-IMEMVCR---KMDDSS---SSQKTEVSTVSSVGP-',
'M-------G----------F-LFNS-V-PWDGPQYHIAPWAFCQTFMGF-VF--GT-PLNAVL-IATYKKLRQPLNYILVNISLGFI-CISVF---VFIS---GYFF-G-VCAEG-FLG-SA-GLVTGWSLAFLAFERFVICKPFGNFRFS-KHSV----VTWIGIGVSI--PPFFGWSRYPEGL-QCSCGPDWYTVGTK---YSE-YYT-FLFICFIPL--LICFSY-QLL-G-ALRAVAAQQQESATTQKAEREVSRMVVMMVGSFCLCYPYAAAMYVNN-R-NHG--LRLVTIPAFFSKSACVYNPIIY-FMNKQFAC-IMEMVCR---KMDDSS---SSQKTEVSTVSSVGP-',
'M-SKMP----EE-------F-LFNSLV-PWDGPQYHLAPWVF-QAFMGF-VF--GT-PLNAVL-VATYRKLRQPLNYILVNVSLGFI-CFSVF---VFISC--GYFF-G-VCAEA-FLG-SAAGLVTGWSLAFLAFERYIICKPFGNFRFS-KHAIA---VTWIGIGVSI--PPFFGWSRFPEGL-QCSCGPDWYTVGTK---YSE-YYT-FLFICYIPL--LICFSY-QLL--RALRAVAAQQQESASTQKAEREVSHMVVVMVGSFCVCYPYAAAMYVNN-R-NHG--LRLVTIPAFFSKSACIYNPIIY-FMNKQFAC-IMEMVC--G-KMDESS---SSQKTEVSTVSSVGP-',
'M-SKM-----E--------F-LFNSLV-PWDGPQYHLAPWAF-QAFMGF-VF--GT-PLNAVL-VATYRKLRQPLNYILVNVSLGFI-CFSVFI--VFISC--GYFF-G-VCAEA-FLGCTA-GLVTGWSLAFLAFERYIICKPFGNFRFS-KHAV----VTWIGIGVSI--PPFFGWSRFPEGL-QCSCGPDWYTVGTK---YSE-YYT-FLFICYIPL--LICFSY-QLL-G-ALRAVAAQQQESASTQKAEREVSHMVVVMVGSFCLCYPYAAAMYVNN-R-NHG--LRLVTIPAFFSKSACVYNPIIY-FMNKQFAC-IMEMVC--G-KMDESS---SSQKTEVSTVSSVGP-',
'------------MDAWAVQFG--NS-V-PFEGEQYHIAPWAF-QAFMGF-VF--GT-PMNGVLFV-TYKKLRQPLNYILVNISLGFID-FSV--SQVFV-CAAGYYFLGTLCAEAA-MG-SA-GLVTGWSLAVLAFERYVICKPFGSFKFQGQ-AV-GAVVTWI-IGTACATPPFFGWSRYPEGLGT-ACGPDWYT---KSEEYSESY-T-FLLICFMPM-III-FSY-QLL-G-ALRAVAAQQAESESTQKAEREVSRMVVVMVGSFVLCYPYAVAMYANSDEPNK--Y-RLVAIPAFFSKSSCVYNPLIY-FMNKQFAC-IMETV--FGKKIDESS-EVSS-KTE--T-SSV--A']
    #print(align.count_differences_in_seq(seq1, seq2))