from constants import *
from numpy import array, delete, zeros, ones, transpose, shape
from numpy import max as np_max
from upgma import Upgma


class Alignement():

    def __init__(self):
        self.upgma = Upgma()
        self.clades_names = []  # Liste des noms de séquences (en seqX)
        self.matrice_distance = []  # Matrice de distance
        self.dir_matrix = []  # Matrice des directions (Après 2e NW)
        self.dict_tree = {}  # Dico après premier UPGMA arbre {'seqX' : distance}
        # Dico avec les séquences pour alignement {'seqX' :  AGTCGAT...}
        self.dict_sequence_tree = {}
        self.alignement1 = []  # Premier alignement
        self.alignement2 = []  # Deuxième alignement si nécessaire
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
        list_seq = []
        list_dist = []
        list_boolean = []
        # On trie le dico avec en clé les séquences et en valeur la distance (les plus proches en premiers)
        for k, v in sorted(self.dict_tree.items(), key=lambda x: x[1]):
            del self.dict_tree[k]
            self.dict_tree[k] = v
        # On place les clés et valeurs dans deux listes séparées (c'est trié) :
        for key, value in self.dict_tree.items():
            list_seq.append(key)
            list_dist.append(value)
        for element in list_dist:
            list_boolean.append(False)
        # On sort les séquences à merge du bon dico qui les contient :

        print(list_dist)
        for i in range(len(list_dist)):
            print(list_boolean)
            if list_boolean[i] == False:
                print("IIIICI : ", i)
                print("self_alignement1 : ", self.alignement1)
                print("\n Etape 0 : ", i)
                # Si c'est la dernière séquence, on l'ajoute à l'alignement:
                if i == len(list_dist)-1:
                    print("coucou")
                    espece_solo = self.dict_sequence_tree[list_seq[i]]
                    self.alignement1 = self.needleman_wunsch_profile(self.alignement1, espece_solo, BLOSUM62)
                    list_boolean[i] = True
                    print("FIN : ", self.alignement1)
                    return self.alignement1
                # Si les deux première distances (triées) sont identiques :
                elif list_dist[i] == list_dist[i+1]:
                    print("\n Etape 0 bis : ", i)
                    if self.alignement1 == []:
                        print("\n Etape 1 : ", i)
                        espece1 = self.dict_sequence_tree[list_seq[i]]
                        espece2 = self.dict_sequence_tree[list_seq[i+1]]
                        self.alignement1 = self.needleman_wunsch_profile(espece1, espece2, BLOSUM62)
                    else:
                        print("\n Etape 1bis : ", i)
                        espece1 = self.dict_sequence_tree[list_seq[i]]
                        espece2 = self.dict_sequence_tree[list_seq[i+1]]
                        self.alignement2 = self.needleman_wunsch_profile(espece1, espece2, BLOSUM62)
                    list_boolean[i] = True
                    list_boolean[i+1] = True
                elif list_dist[i] != list_dist[i+1]:
                    print("\n Etape 0ter : ", i)
                    self.alignement2 = [self.dict_sequence_tree[list_seq[i]]]
                    list_boolean[i] = True
                if self.alignement1 != [] and self.alignement2 != []:
                    print("\n Etape X : ", i)
                    self.alignement1 = self.needleman_wunsch_profile(self.alignement1, self.alignement2, BLOSUM62)
                    
        #         i += 2
        #         print("\n Etape 1 : ", i)
        #         if i in range(len(list_dist)-2):
        #             if list_dist[i] == list_dist[i+1]:
        #                 espece3 = self.dict_sequence_tree[list_seq[i]]
        #                 espece4 = self.dict_sequence_tree[list_seq[i+1]]
        #                 i += 2
        #                 print("\n Etape 2 : ", i)
        #                 self.alignement2 = self.needleman_wunsch_profile(
        #                     espece3, espece4, BLOSUM62)
        #             if self.alignement2 != []:
        #                 self.alignement1 = self.needleman_wunsch_profile(
        #                     self.alignement1, self.alignement2, BLOSUM62)
        #                 self.alignement2 = []
        #     elif list_dist[i] != list_dist[i+1]:
        #         print("\n Etape 3 : ", i)
        #         espece_solo = self.dict_sequence_tree[list_seq[i]]
        #         if self.alignement1 == []:
        #             if list_dist[i+1] == list_dist[i+2]:
        #                 espece1 = self.dict_sequence_tree[list_seq[i+1]]
        #                 espece2 = self.dict_sequence_tree[list_seq[i+2]]
        #                 self.alignement1 = self.needleman_wunsch_profile(espece1, espece2, BLOSUM62)
        #                 self.alignement1 = self.needleman_wunsch_profile(self.alignement1, espece_solo, BLOSUM62)
        #                 i += 2
        #         else:
        #             self.alignement1 = self.needleman_wunsch_profile(
        #                     self.alignement1, espece_solo, BLOSUM62)
        print("apres")
        return self.alignement1




if __name__ == "__main__":
    align = Alignement()
    print(align.multiple_alignement())
